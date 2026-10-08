// FEBioAlignerPairs: the aligner-crown contact of ofj4d/solve4d.py as an FEBio 4 nonlinear constraint
// ("aligner-pairs"), for the FEBio backend ofj4d/febio4d.py (docs/FINDINGS.md section 31).
//
// FEBio's own contact re-projects every contact point along its normal at every iteration; on this
// appliance the partner facet jumps (embrasures, curved crowns) and the resulting force jumps invert
// elements under every setting tried. solve4d instead keeps the pairing of each aligner node fixed for
// the whole solve (small sliding). This constraint implements exactly that law, so that FEBio differs
// from solve4d only in the continuum (finite deformation, neo-Hookean aligner, finite-strain PDL):
//
//   pairs: slave node s -> master triangle (m0, m1, m2) with barycentric weights w, the master normal n
//          at the start of the solve, and the reference gap gref (solve4d: gap_now - gap_pre + clearance)
//   rel  = u_s - sum_i w_i u_mi                       (u = displacement since the start of the solve)
//   cn   = ramp(t) * gref + n . rel                   (ramp on a load curve: 0 = exactly passive)
//   lam  = K_N max(-cn, 0)
//   t    = rel - (n . rel) n,  k_sec = K_T if K_T |t| <= mu lam, else mu lam / |t|   (secant Coulomb)
//   force on the slave node  f = lam n - k_sec t,  on master node i  -w_i f
//   tangent  B^T (K_N n n + k_sec (I - n n)) B  for cn <= 0 (so a passive appliance is seated), and the
//            normal part K_N n n for pairs separated by less than stab_band (no force): as the interference
//            ramps in, pairs separate and a region of the aligner can lose its seating, leaving a
//            near-singular tangent; the band keeps it seated without changing the equilibrium
//
// Nodes attached to rigid bodies are handled by FEBio's residual / linear-system assembly (node lists).
#include <FECore/sdk.h>
#include <FECore/FENLConstraint.h>
#include <FECore/FELinearSystem.h>
#include <FECore/FEGlobalMatrix.h>
#include <FECore/FEModel.h>
#include <FECore/FEMesh.h>
#include <FECore/log.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

class FEAlignerPairs : public FENLConstraint
{
public:
	explicit FEAlignerPairs(FEModel* fem) : FENLConstraint(fem) {}

	std::string m_file;
	double m_kn = 2000.0, m_kt = 2000.0, m_mu = 0.2, m_ramp = 1.0, m_band = 0.0;

	struct Pair { int s; int m[3]; double w[3]; vec3d n; double gref; };
	std::vector<Pair> m_pairs;
	int m_dof[6] = { -1, -1, -1, -1, -1, -1 };

	bool Init() override
	{
		FEModel& fem = *GetFEModel();
		const char* names[6] = { "x", "y", "z", "Ru", "Rv", "Rw" };
		for (int k = 0; k < 6; ++k) m_dof[k] = fem.GetDOFIndex(names[k]);
		FEMesh& mesh = GetMesh();
		std::ifstream in(m_file);
		if (!in) { feLogError("aligner-pairs: cannot open %s", m_file.c_str()); return false; }
		m_pairs.clear();
		std::string line;
		while (std::getline(in, line))
		{
			if (line.empty() || line[0] == '#') continue;
			std::replace(line.begin(), line.end(), ',', ' ');
			std::istringstream ss(line);
			double id[4], v[7];
			for (double& x : id) ss >> x;
			for (double& x : v) ss >> x;
			if (!ss) { feLogError("aligner-pairs: bad line in %s", m_file.c_str()); return false; }
			Pair p;
			p.s = mesh.FindNodeIndexFromID((int)id[0]);
			for (int k = 0; k < 3; ++k) p.m[k] = mesh.FindNodeIndexFromID((int)id[k + 1]);
			for (int k = 0; k < 3; ++k) p.w[k] = v[k];
			p.n = vec3d(v[3], v[4], v[5]); p.n.unit();
			p.gref = v[6];
			if (p.s < 0 || p.m[0] < 0 || p.m[1] < 0 || p.m[2] < 0) { feLogError("aligner-pairs: unknown node id"); return false; }
			m_pairs.push_back(p);
		}
		feLog("aligner-pairs: %d pairs from %s, K_N %g, K_T %g, mu %g, stab_band %g\n", (int)m_pairs.size(), m_file.c_str(), m_kn, m_kt, m_mu, m_band);
		return FENLConstraint::Init();
	}

	// kinematics and forces of one pair
	void Eval(const Pair& p, double& cn, double& lam, double& ksec, vec3d& t) const
	{
		FEMesh& mesh = const_cast<FEAlignerPairs*>(this)->GetMesh();
		auto u = [&](int i) { const FENode& nd = mesh.Node(i); return nd.m_rt - nd.m_r0; };
		vec3d rel = u(p.s) - (u(p.m[0]) * p.w[0] + u(p.m[1]) * p.w[1] + u(p.m[2]) * p.w[2]);
		double un = p.n * rel;
		cn = m_ramp * p.gref + un;
		lam = (cn < 0 ? -m_kn * cn : 0.0);
		t = rel - p.n * un;
		double tm = t.norm();
		ksec = 0.0;
		if (cn <= 0)
		{
			ksec = (m_kt * tm <= m_mu * lam ? m_kt : m_mu * lam / std::max(tm, 1e-12));
			ksec = std::min(std::max(ksec, 0.0), m_kt);
		}
	}

	void Layout(const Pair& p, std::vector<int>& en, std::vector<int>& lm) const
	{
		FEMesh& mesh = const_cast<FEAlignerPairs*>(this)->GetMesh();
		en = { p.s, p.m[0], p.m[1], p.m[2] };
		lm.resize(12);
		for (int a = 0; a < 4; ++a)
		{
			const std::vector<int>& id = mesh.Node(en[a]).m_ID;
			for (int k = 0; k < 3; ++k) lm[3 * a + k] = id[m_dof[k]];
		}
	}

	void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override
	{
		std::vector<int> en, lm; std::vector<double> fe(12);
		double Ftot = 0; int nact = 0; vec3d Fsum(0, 0, 0);
		for (const Pair& p : m_pairs)
		{
			double cn, lam, ksec; vec3d t;
			Eval(p, cn, lam, ksec, t);
			if (cn > 0) continue;
			vec3d f = p.n * lam - t * ksec;
			double c[4] = { 1.0, -p.w[0], -p.w[1], -p.w[2] };
			for (int a = 0; a < 4; ++a) { fe[3 * a] = c[a] * f.x; fe[3 * a + 1] = c[a] * f.y; fe[3 * a + 2] = c[a] * f.z; }
			Layout(p, en, lm);
			R.Assemble(en, lm, fe);
			if (lam > 0) { ++nact; Fsum += p.n * lam; }
		}
		m_nact = nact; m_Fn = Fsum.norm();
	}

	void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override
	{
		std::vector<int> en, lm;
		for (const Pair& p : m_pairs)
		{
			double cn, lam, ksec; vec3d t;
			Eval(p, cn, lam, ksec, t);
			if (cn > m_band) continue;
			if (cn > 0) ksec = 0.0;                         // stabilisation only: normal stiffness, no friction
			double k[3][3];
			const double nv[3] = { p.n.x, p.n.y, p.n.z };
			for (int i = 0; i < 3; ++i)
				for (int j = 0; j < 3; ++j)
					k[i][j] = m_kn * nv[i] * nv[j] + ksec * ((i == j ? 1.0 : 0.0) - nv[i] * nv[j]);
			double c[4] = { 1.0, -p.w[0], -p.w[1], -p.w[2] };
			FEElementMatrix ke; ke.resize(12, 12); ke.zero();
			for (int a = 0; a < 4; ++a)
				for (int b = 0; b < 4; ++b)
					for (int i = 0; i < 3; ++i)
						for (int j = 0; j < 3; ++j)
							ke[3 * a + i][3 * b + j] = c[a] * c[b] * k[i][j];
			Layout(p, en, lm);
			ke.SetNodes(en);
			ke.SetIndices(lm);
			LS.Assemble(ke);
		}
	}

	void BuildMatrixProfile(FEGlobalMatrix& M) override
	{
		FEMesh& mesh = GetMesh();
		std::vector<int> lm(24);
		for (const Pair& p : m_pairs)
		{
			int nodes[4] = { p.s, p.m[0], p.m[1], p.m[2] };
			for (int a = 0; a < 4; ++a)
			{
				const std::vector<int>& id = mesh.Node(nodes[a]).m_ID;
				for (int k = 0; k < 6; ++k) lm[6 * a + k] = (m_dof[k] >= 0 ? id[m_dof[k]] : -1);
			}
			M.build_add(lm);
		}
	}

	bool Augment(int naug, const FETimeInfo& tp) override
	{
		feLog("aligner-pairs: %d pairs in compression, |sum lam n| = %g N (ramp %g)\n", m_nact, m_Fn, m_ramp);
		return true;
	}

	int m_nact = 0; double m_Fn = 0;

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(FEAlignerPairs, FENLConstraint)
	ADD_PARAMETER(m_file, "pairs_file");
	ADD_PARAMETER(m_kn, "penalty");
	ADD_PARAMETER(m_kt, "tangential_penalty");
	ADD_PARAMETER(m_mu, "fric_coeff");
	ADD_PARAMETER(m_ramp, "ramp")->SetFlags(FE_PARAM_ADDLC | FE_PARAM_VOLATILE);
	ADD_PARAMETER(m_band, "stab_band");
END_FECORE_CLASS();

FECORE_PLUGIN int GetSDKVersion() { return FE_SDK_VERSION; }

FECORE_PLUGIN void PluginInitialize(FECoreKernel& febio)
{
	FECoreKernel::SetInstance(&febio);
	REGISTER_FECORE_CLASS(FEAlignerPairs, "aligner-pairs");
}
