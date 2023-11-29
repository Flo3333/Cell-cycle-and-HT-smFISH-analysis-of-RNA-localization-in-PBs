# -*- coding: utf-8 -*-
# Author: Floric Slimani <floric.slimani@live.fr>


"""
This subpackage handles data integrity amongst our custom pandas DataFrames. 
"""

from .datashape import check_id, check_samedatashape, check_expectedcolumns, check_isnotempty
from .datashape import is_empty, is_primarykey, is_contained, is_primary
from .datashape import has_samedatashape, has_column, has_id
from .Errors import MissingColumnError, MergeError