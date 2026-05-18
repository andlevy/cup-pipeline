class BackendError(Exception):
    """
        If the system is not in a valid state to run the code
        e.g. reference files are not accessible, errors with alpha missense or ensembl vep
        Use if admin needs to debug the backend
    """

    def __init__(self, message="Backend setup incorrectly"):
        super().__init__(message)
        self.error_type = "BackendError"


class FrontendError(Exception):
    """
        If there is something unexpected about the input that the user didn't have control over
        e.g. the files passed had incorrect permissions
        Use if admin needs to debug the frontend
    """

    def __init__(self, message="Frontend/cli use setup incorrectly"):
        super().__init__(message)
        self.error_type = "FrontendError"


class IncorrectInputError(Exception):
    """
        If the user input is incorrect
        e.g. the files are missing or in the wrong format
        Use if the user needs to change their input
        This error message must be suitable for website display
    """

    def __init__(self, message="Input files are incorrect"):
        super().__init__(message)
        self.error_type = "IncorrectInputError"
