import feather
import pandas as pd

# read feather file with rownames stored in 1st column to a pandas df with rownames
# properly stored as index
def import_feather(fn):
    try:
        df = feather.read_dataframe(fn)
    except OSError as err:
        print("OS error: {0}".format(err))
        
    df.set_index(df.iloc[:,0], inplace=True)
    df.index.name=None
    df.drop(df.columns[0],axis=1,inplace=True)
    return df 
    
# feather drops index/rownames for some reason, so save them as their own column
def export_feather(df, fn):
    df.insert(loc=0, column="rowname", value=df.index)
    try:
        feather.write_dataframe(df, fn)
    except OSError as err:
        print("OS error: {0}".format(err))
