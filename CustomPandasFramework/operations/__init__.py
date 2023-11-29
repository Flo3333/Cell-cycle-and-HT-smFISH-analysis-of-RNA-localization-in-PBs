# -*- coding: utf-8 -*-
# Author: Floric Slimani <floric.slimani@live.fr>


"""
This subpackage handles data operations amongst our custom pandas DataFrames. 
"""

from .frame_operations import add_data, foreign_key, get_referencement_relation, set_missingcolumns_toNA, get_missing_column, keep_columns