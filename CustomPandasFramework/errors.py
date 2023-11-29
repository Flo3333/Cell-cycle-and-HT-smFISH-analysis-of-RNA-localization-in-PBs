"""
This sub module contains dataframes integrity related errors
"""

class DataframeIntegrityError(Exception):
    """
    General Exception class for DataFrame integrity erros.
    """
    pass


class EmptyFrameError(Exception):
    pass



### Primary Key Errors

class PrimaryKeyError(DataframeIntegrityError):
    """
    Exception class for non valid id issues.
    """
    pass

class NoIdError(PrimaryKeyError):
    pass

class IdIsNotPrimaryError(PrimaryKeyError):
    pass

class MissingIdsError(PrimaryKeyError):
    pass


### Data Shape Errors
class DataShapeError(DataframeIntegrityError):
    """
    Execption class raised for issues with frames shapes.
    """
    pass

class MissingColumnsError(DataShapeError):
    pass

class MissingIndexError(DataShapeError):
    pass

class DataLengthError(DataShapeError):
    pass

