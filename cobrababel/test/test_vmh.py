import pytest
import cobrababel


class TestVirtualMetabolicHuman:
    def test_create_model(self):
        model = cobrababel.create_cobra_model_from_agora_model('Bacteroides_thetaiotaomicron_VPI_5482')
        assert len(model.metabolites) == 1146
        assert len(model.reactions) == 1362
        assert len(model.genes) == 840
        solution = model.optimize()
        assert solution.f == pytest.approx(0.394520)

    def test_create_recon2(self):
        recon2 = cobrababel.create_cobra_model_from_vmh_recon2()
        assert len(recon2.metabolites) >= 5063
        assert len(recon2.reactions) >= 7440
        assert len(recon2.genes) >= 2140
        solution = recon2.optimize()
        assert solution.f == pytest.approx(3.198056)

