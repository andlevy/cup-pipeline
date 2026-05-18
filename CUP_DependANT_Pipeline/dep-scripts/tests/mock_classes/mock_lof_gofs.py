from app.lofs_and_gofs.ensembl_vep import EnsemblVep


class MockEnsemblVep(EnsemblVep):

    def run_ensembl_vep(self):
        return self.data.vep_results
