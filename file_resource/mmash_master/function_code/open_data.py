def create_dataset(users,file_name,path,replace_na=True):
    
    """
    Returns MMASH datafarame.
    Parameters
    ---------
    nn_intervals : list
        List of Normal to Normal Interval.
    Returns
    ---------
    nni_tmstp : list
        list of time intervals between first NN-interval and final NN-interval.
    """
    
    import pandas
    import numpy
    
    df_concat = pandas.DataFrame()
    for user in users:
        try:
            df = pandas.read_csv('%s/%s/%s.csv'%(path,user,file_name))
            df['user'] = user
            df_concat = pandas.concat([df_concat,df])
        except:
            print('NO data for %s'%user)
            pass
    
    del df_concat['Unnamed: 0']
    df_concat = df_concat.set_index('user')

    if replace_na == True:
        df_concat = df_concat.replace(0,numpy.nan)

    return(df_concat)