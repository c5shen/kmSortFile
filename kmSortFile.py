import sys
import numpy as np
import pandas as pd
from scipy import stats
from cuzcatlan import elemental

def kmSortFile(num_clusters, name):

    P_VALUE = 0.005

    # input a number num_clusters to indicate number of clusters assigned
    # input a name for the basis of cluster files name
    #num_clusters = int(sys.argv[1])
    #name = sys.argv[2]

    # original data with combined clusters
    original = open(name+'.gct')

    clusters = []
    for i in range(1,num_clusters+1):
        clusters.append(open(name + '-' + str(i) + '.gct'))

    temp = original.read().split('\n')
    header = temp[0]
    info = temp[1]
    num_rows = int(temp[1].split('\t')[0])
    num_columns = int(temp[1].split('\t')[1])

    # store the samples of the cluster, and also column labels for pandas
    # dataframe construction
    column_labels = temp[2].split('\t')
    sample_names = temp[2].split('\t')[2:]

    temp_data = []      # temporarily holds data for future numpy format
    description = []    # stores descriptions for each probe
    probe_names = []    # stores names of each probe

    # pull down the data to different chunks.
    # 1) description
    # 2) probe names
    # 3) data part -> matrix
    for i in range(3, len(temp)):
        if temp[i] != '':
            processing = temp[i].split('\t')
            temp_data.append(list(processing[2:]))
            description.append(processing[1])
            probe_names.append(processing[0])
            del processing
    matrix = np.array(temp_data, dtype=float)

    # get clusters info (num of rows, starting location in matrix)
    clusters_info = []
    starting_loc = 0
    for i in range(num_clusters):
        temp_num = int(clusters[i].read().split('\n')[1].split('\t')[1])
        clusters_info.append((starting_loc, temp_num))
        starting_loc += temp_num

    # use a dictionary to store all uniquely expressed clusters
    # the key is index of probe, and the value is the index of cluster
    # a tuple (_which_cluster_, _p_value_) is stored in each entry
    diff_expression = {}

    # Find uniquely expressed/repressed genes of the provided clusters
    for i in range(num_clusters):
        start = clusters_info[i][0]
        length = clusters_info[i][1]

        # for each gene, flag the current cluster as test, and others as
        # compared calculate the mean for each and do a two-sample t-test to
        # see if the cluster gene expression is uniquely expressed
        for j in range(num_rows):
            test = np.copy(matrix[j][start:start+length])
            compared = np.copy(matrix[j])
            for item in test:
                compared = np.delete(compared, item)
            
            # if the expression level is statistically significant (0.05),
            # both tails, then include this in the different expressed clusters
            p_value = stats.ttest_ind(test, compared)[1]
            if (p_value < P_VALUE):

                # if there is already another cluster has the same differently
                # expressed gene, check the p_value and the smaller one wins
                if j in diff_expression:
                    old_p_value = diff_expression[j][1]
                    if p_value < old_p_value:
                        diff_expression[j] = (i, p_value)
                
                # else, we annotate this diff expressed gene to the
                # corresponding cluster
                else:
                    diff_expression[j] = (i, p_value)

    # output to a file named after the original name, with 'sorted' tagged
    # behind
#    output = open(name + '-sorted.gct', 'w')

    # write the header and # of rows/columns, also write the names of sample
#    output.write(header + '\n' + info + '\n')
#    output.write("Name" + '\t' + "Description" + '\t' + 
#                      '\t'.join(sample_names) + '\n')


    # construct a int array containing length of num of rows, start writing in
    # rows that are differently expressed by each cluster one by one,
    # simultaneously remove them from this array. Finally, write the rest rows
    # IN ORDER to the file
    rest_rows = [x for x in range(num_rows)]

    # iterate through the whole dictionary for each cluster, so we can print
    # out the corresponding differently expressed genes in cluster order
#    for k in range(num_clusters):
#        for key,value in diff_expression.items():
#            if value[0] == k:
#                output.write(probe_names[key] + '\t' + description[key] + '\t')
#                for i in range(len(matrix[key])):
#                    output.write(str(matrix[key][i]) + '\t')
#                output.write('\n')
#                rest_rows.remove(key)   # remove the added ones from rest_rows

#    for i in range(len(rest_rows)):
#        output.write(probe_names[rest_rows[i]] + '\t' + 
#                         description[rest_rows[i]] + '\t')
#        for j in range(len(matrix[rest_rows[i]])):
#            output.write(str(matrix[rest_rows[i]][j]) + '\t')
#        output.write('\n')
        
#    output.close()

    df_to_process = pd.read_table(name+'.gct', skiprows=2)
    new_order = []
    # rearrange the rows to pass the matrix for an np ndarray, then construct
    # pandas DataFrame
    for key, value in sorted(diff_expression.items(), key=lambda item :(item[1][0],item[0]) ):
        #print("%s: %s" % (key,value))
        rest_rows.remove(key)
        new_order.append(key)
    for item in rest_rows:
        new_order.append(item)

    #print(df_to_process.index.values)
    df_to_process = df_to_process.reindex(new_order)
    #print(df_to_process)

    elemental.df2gct(df_to_process, 2, True, name+'-sorted.gct', False)

kmSortFile(int(sys.argv[1]), sys.argv[2])
