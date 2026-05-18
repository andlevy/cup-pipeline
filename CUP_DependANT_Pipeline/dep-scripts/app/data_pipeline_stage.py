class DataPipelineStage:

    def __init__(self, data):
        self.data = data
        return

    def execute(self):
        raise NotImplementedError("Execute method must be overridden in subclasses")
