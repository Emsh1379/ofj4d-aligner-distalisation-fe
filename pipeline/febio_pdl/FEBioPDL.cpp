// FEBioPDL: the bilinear PDL of Mao et al. (Vollmer's "double linear" curve) as an FEBio 4
// material plug-in, for the FEBio backend of the staging model (ofj4d/febio4d.py).
//
// Same law as ofj4d/solve4d.py -- a secant modulus driven by the von Mises equivalent
// strain,  E_s = E1  (e <= ek),  (E1 ek + E2 (e - ek)) / e  (e > ek)  -- but on
// finite-strain kinematics: e is the equivalent Green-Lagrange strain and the stress is
// compressible neo-Hookean with modulus E_s. In the small-strain limit it reduces exactly
// to the custom solver's material. The tangent is the neo-Hookean tangent at the INCREMENTAL
// modulus (E1 below the knee, E2 above), not at E_s: the secant modulus under-estimates the stiffness
// above the knee, and Newton then overshoots and inverts soft PDL elements (docs/FINDINGS.md
// section 31). The stress is unchanged, so the converged solution is too.
#include <FECore/sdk.h>
#include <FEBioMech/FEElasticMaterial.h>
#include <cmath>

class FEPDLBilinear : public FEElasticMaterial
{
public:
	explicit FEPDLBilinear(FEModel* fem) : FEElasticMaterial(fem) {}

	double m_E1 = 0.05, m_E2 = 0.22, m_ek = 0.075, m_v = 0.45;

	double Esec(FEElasticMaterialPoint& pt) const
	{
		mat3ds E = pt.Strain();                              // Green-Lagrange
		double tr = E.tr()/3.0;
		mat3ds d = E - mat3dd(tr);
		double e = std::sqrt(std::max(2.0/3.0*d.dotdot(d), 0.0));
		if (e <= m_ek) return m_E1;
		return (m_E1*m_ek + m_E2*(e - m_ek))/e;
	}

	mat3ds Stress(FEMaterialPoint& mp) override
	{
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		double J = pt.m_J, E = Esec(pt);
		double lam = m_v*E/((1 + m_v)*(1 - 2*m_v)), mu = 0.5*E/(1 + m_v);
		mat3ds b = pt.LeftCauchyGreen();
		mat3dd I(1);
		return (b - I)*(mu/J) + I*(lam*std::log(J)/J);
	}

	// incremental modulus of the bilinear curve at this strain
	double Etan(FEElasticMaterialPoint& pt) const
	{
		mat3ds E = pt.Strain();
		mat3ds d = E - mat3dd(E.tr()/3.0);
		double e = std::sqrt(std::max(2.0/3.0*d.dotdot(d), 0.0));
		return (e <= m_ek ? m_E1 : m_E2);
	}

	tens4ds Tangent(FEMaterialPoint& mp) override
	{
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		double J = pt.m_J, E = Etan(pt);
		double lam = m_v*E/((1 + m_v)*(1 - 2*m_v)), mu = 0.5*E/(1 + m_v);
		mat3dd I(1);
		return dyad1s(I)*(lam/J) + dyad4s(I)*(2*(mu - lam*std::log(J))/J);
	}

	double StrainEnergyDensity(FEMaterialPoint& mp) override
	{
		FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();
		double J = pt.m_J, lnJ = std::log(J), E = Esec(pt);
		double lam = m_v*E/((1 + m_v)*(1 - 2*m_v)), mu = 0.5*E/(1 + m_v);
		return mu*((pt.LeftCauchyGreen().tr() - 3)/2.0 - lnJ) + lam*lnJ*lnJ/2.0;
	}

	DECLARE_FECORE_CLASS();
};

BEGIN_FECORE_CLASS(FEPDLBilinear, FEElasticMaterial)
	ADD_PARAMETER(m_E1, "E1");
	ADD_PARAMETER(m_E2, "E2");
	ADD_PARAMETER(m_ek, "eps_knee");
	ADD_PARAMETER(m_v , "v");
END_FECORE_CLASS();

FECORE_PLUGIN int GetSDKVersion() { return FE_SDK_VERSION; }

FECORE_PLUGIN void PluginInitialize(FECoreKernel& febio)
{
	FECoreKernel::SetInstance(&febio);
	REGISTER_FECORE_CLASS(FEPDLBilinear, "pdl bilinear");
}
