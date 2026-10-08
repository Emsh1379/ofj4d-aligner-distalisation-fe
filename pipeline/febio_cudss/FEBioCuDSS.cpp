// FEBioCuDSS: NVIDIA cuDSS sparse direct solver as an FEBio 4 linear-solver plug-in.
//
// Compiled and verified against PARDISO on 2026-09-13 (FEBio 4.13.0, cuDSS 0.8.0):
// identical iterations, node displacements equal within 1e-12 mm, solver 5.1x faster
// on a T4 (README.md). The linear-solver call counters are not updated yet.
//
// FEBio hands the solver a 1-based CSR matrix: CompactSymmMatrix(1) holds the upper
// triangle of a symmetric matrix, CRSSparseMatrix(1) a full unsymmetric one. The
// three FEBio phases map one-to-one onto cuDSS:
//   PreProcess (sparsity known) -> allocate device arrays, create cuDSS objects
//   Factor                      -> copy values, ANALYSIS (first time), FACTORIZATION
//   BackSolve                   -> SOLVE
#include <FECore/sdk.h>
#include <FECore/LinearSolver.h>
#include <FECore/CompactSymmMatrix.h>
#include <FECore/CompactUnSymmMatrix.h>
#include <cuda_runtime_api.h>
#include <cudss.h>
#include <FECore/log.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>

class CuDSSSolver : public LinearSolver
{
public:
	explicit CuDSSSolver(FEModel* fem) : LinearSolver(fem) {}
	~CuDSSSolver() override { Destroy(); }

	SparseMatrix* CreateSparseMatrix(Matrix_Type ntype) override
	{
		switch (ntype)
		{
		case REAL_SYMMETRIC     : m_sym = true;  m_pA = new CompactSymmMatrix(1); break;
		case REAL_UNSYMMETRIC   :
		case REAL_SYMM_STRUCTURE: m_sym = false; m_pA = new CRSSparseMatrix(1); break;
		default: m_pA = nullptr;
		}
		return m_pA;
	}

	bool SetSparseMatrix(SparseMatrix* pA) override
	{
		m_pA = dynamic_cast<CompactMatrix*>(pA);
		m_sym = (dynamic_cast<CompactSymmMatrix*>(pA) != nullptr);
		return (m_pA != nullptr);
	}

	bool PreProcess() override
	{
		if (m_pA == nullptr) return false;
		Release();                                   // sparsity may have changed (contact)
		m_n = m_pA->Rows();
		m_nnz = (int)m_pA->NonZeroes();
		if (m_n == 0) return true;

		if (!ok(cudaMalloc(&d_ptr, sizeof(int)*(m_n + 1))) ||
		    !ok(cudaMalloc(&d_ind, sizeof(int)*m_nnz)) ||
		    !ok(cudaMalloc(&d_val, sizeof(double)*m_nnz)) ||
		    !ok(cudaMalloc(&d_x, sizeof(double)*m_n)) ||
		    !ok(cudaMalloc(&d_b, sizeof(double)*m_n))) return false;
		if (!ok(cudaMemcpy(d_ptr, m_pA->Pointers(), sizeof(int)*(m_n + 1), cudaMemcpyHostToDevice)) ||
		    !ok(cudaMemcpy(d_ind, m_pA->Indices(), sizeof(int)*m_nnz, cudaMemcpyHostToDevice))) return false;

		if (cudssCreate(&m_h) != CUDSS_STATUS_SUCCESS ||
		    cudssConfigCreate(&m_cfg) != CUDSS_STATUS_SUCCESS ||
		    cudssDataCreate(m_h, &m_data) != CUDSS_STATUS_SUCCESS) return false;

		// robustness for symmetric INDEFINITE systems (rigid bodies, ties, penalty contact):
		// iterative refinement, optional diagonal pivoting (FEBIO_CUDSS_IR_STEPS, FEBIO_CUDSS_PIVOT)
		{
			const char* e = std::getenv("FEBIO_CUDSS_IR_STEPS");
			int ir = e ? std::atoi(e) : 2;
			if (ir > 0) cudssConfigSet(m_cfg, CUDSS_CONFIG_IR_N_STEPS, &ir, sizeof(ir));
			const char* p = std::getenv("FEBIO_CUDSS_PIVOT");
			if (p && std::strcmp(p, "diagonal") == 0) { cudssPivotType_t pt = CUDSS_PIVOT_DIAGONAL; cudssConfigSet(m_cfg, CUDSS_CONFIG_PIVOT_TYPE, &pt, sizeof(pt)); }
			if (p && std::strcmp(p, "none") == 0) { cudssPivotType_t pt = CUDSS_PIVOT_NONE; cudssConfigSet(m_cfg, CUDSS_CONFIG_PIVOT_TYPE, &pt, sizeof(pt)); }
			const char* c = std::getenv("FEBIO_CUDSS_RESCHECK");
			m_rescheck = c ? std::atof(c) : 1e-8;
		}

		const cudssIndexBase_t base = (m_pA->Offset() == 1 ? CUDSS_BASE_ONE : CUDSS_BASE_ZERO);
		if (cudssMatrixCreateCsr(&m_A, m_n, m_n, m_nnz, d_ptr, nullptr, d_ind, d_val,
		        CUDSS_R_32I, CUDSS_R_32I, CUDSS_R_64F,       // offsets, indices, values
		        m_sym ? CUDSS_MTYPE_SYMMETRIC : CUDSS_MTYPE_GENERAL,
		        m_sym ? CUDSS_MVIEW_UPPER : CUDSS_MVIEW_FULL, base) != CUDSS_STATUS_SUCCESS) return false;
		if (cudssMatrixCreateDn(&m_x, m_n, 1, m_n, d_x, CUDSS_R_64F, CUDSS_LAYOUT_COL_MAJOR) != CUDSS_STATUS_SUCCESS ||
		    cudssMatrixCreateDn(&m_b, m_n, 1, m_n, d_b, CUDSS_R_64F, CUDSS_LAYOUT_COL_MAJOR) != CUDSS_STATUS_SUCCESS)
			return false;
		m_analysed = false;
		return true;
	}

	bool Factor() override
	{
		if (m_n == 0) return true;
		if (!ok(cudaMemcpy(d_val, m_pA->Values(), sizeof(double)*m_nnz, cudaMemcpyHostToDevice))) return false;
		if (!m_analysed)
		{
			if (cudssExecute(m_h, CUDSS_PHASE_ANALYSIS, m_cfg, m_data, m_A, m_x, m_b) != CUDSS_STATUS_SUCCESS)
				return false;
			m_analysed = true;
		}
		if (cudssExecute(m_h, CUDSS_PHASE_FACTORIZATION, m_cfg, m_data, m_A, m_x, m_b) != CUDSS_STATUS_SUCCESS) return false;
		// diagnostics: info != 0 flags a (near-)zero pivot; npivots counts perturbed pivots
		int info = 0, npiv = 0; size_t w = 0;
		cudssDataGet(m_h, m_data, CUDSS_DATA_INFO, &info, sizeof(info), &w);
		cudssDataGet(m_h, m_data, CUDSS_DATA_NPIVOTS, &npiv, sizeof(npiv), &w);
		++m_nfact;
		if (info != 0 || npiv != 0) feLog("cudss: factorisation %d: info %d, perturbed pivots %d (n = %d)\n", m_nfact, info, npiv, m_n);
		return true;
	}

	bool BackSolve(double* x, double* y) override
	{
		if (m_n == 0) return true;
		if (!ok(cudaMemcpy(d_b, y, sizeof(double)*m_n, cudaMemcpyHostToDevice))) return false;
		if (cudssExecute(m_h, CUDSS_PHASE_SOLVE, m_cfg, m_data, m_A, m_x, m_b) != CUDSS_STATUS_SUCCESS)
			return false;
		if (!ok(cudaMemcpy(x, d_x, sizeof(double)*m_n, cudaMemcpyDeviceToHost))) return false;   // synchronises
		// true relative residual on the host, with FEBio's own matrix product
		if (m_rescheck > 0)
		{
			std::vector<double> r(m_n);
			m_pA->mult_vector(x, r.data());
			double nb = 0, nr = 0;
			for (int i = 0; i < m_n; ++i) { nb += y[i]*y[i]; nr += (y[i] - r[i])*(y[i] - r[i]); }
			double rel = (nb > 0 ? std::sqrt(nr/nb) : std::sqrt(nr));
			if (!std::isfinite(rel) || rel > m_rescheck)
				feLog("cudss: solve after factorisation %d: relative residual %.3e (n = %d)\n", m_nfact, rel, m_n);
		}
		return true;
	}

	void Destroy() override { Release(); }

private:
	static bool ok(cudaError_t e) { return e == cudaSuccess; }
	void Release()
	{
		if (m_A) cudssMatrixDestroy(m_A);
		if (m_x) cudssMatrixDestroy(m_x);
		if (m_b) cudssMatrixDestroy(m_b);
		if (m_data) cudssDataDestroy(m_h, m_data);
		if (m_cfg) cudssConfigDestroy(m_cfg);
		if (m_h) cudssDestroy(m_h);
		for (void* p : {(void*)d_ptr, (void*)d_ind, (void*)d_val, (void*)d_x, (void*)d_b}) if (p) cudaFree(p);
		m_A = m_x = m_b = nullptr; m_data = nullptr; m_cfg = nullptr; m_h = nullptr;
		d_ptr = d_ind = nullptr; d_val = d_x = d_b = nullptr;
		m_analysed = false;
	}

	CompactMatrix* m_pA = nullptr;
	bool m_sym = true, m_analysed = false;
	double m_rescheck = 1e-8;
	int m_nfact = 0;
	int m_n = 0, m_nnz = 0;
	cudssHandle_t m_h = nullptr;
	cudssConfig_t m_cfg = nullptr;
	cudssData_t m_data = nullptr;
	cudssMatrix_t m_A = nullptr, m_x = nullptr, m_b = nullptr;
	int* d_ptr = nullptr; int* d_ind = nullptr;
	double* d_val = nullptr; double* d_x = nullptr; double* d_b = nullptr;
};

FECORE_PLUGIN int GetSDKVersion() { return FE_SDK_VERSION; }

FECORE_PLUGIN void PluginInitialize(FECoreKernel& febio)
{
	FECoreKernel::SetInstance(&febio);
	REGISTER_FECORE_CLASS(CuDSSSolver, "cudss");
}
