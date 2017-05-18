import cobrababel


class TestMetaNetX:
    def test_create_universal(self):
        universal = cobrababel.create_metanetx_universal_model()
        assert universal.id == 'metanetx_universal'
        assert len(universal.reactions) >= 31283
        assert len(universal.metabolites) >= 132812
