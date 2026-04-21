import os
import json
import pandas as pd
import argparse

def merge_dicts(dict1, dict2):
    merged_dict = {}
    n = 0 
    for key in set(dict1.keys()).union(dict2.keys()):
        ##print(n)
        ##print(key)
        if key in dict1 and key in dict2:
            tmp_dict = {}
            list_of_dicts = [dict1.get(key), dict2.get(key)]
            #print(list_of_dicts)
            for single_dict in list_of_dicts:
                #print(single_dict)
                for k, v in single_dict.items():
                    if k in tmp_dict:
                        tmp_dict[k] = int(tmp_dict[k]) + int(v)
                    else:
                        tmp_dict[k] = v
            
            merged_dict[key] = tmp_dict
            #print(merged_dict)
            #break
        elif key in dict1:
            merged_dict[key] = dict1[key]
        else:
            merged_dict[key] = dict2[key]
        ##print(merged_dict)
        ##print("--")
        ##n = n+1
        ##if n>10:
        ##    break
    return merged_dict

# check
#dict1 = {"t1": {"t11":1},"t2":{'t21': 1}}
#dict2 = {"t2":{'t21': 1, "a":1},"t3":{"t31":3}}
#merged_dict = merge_dicts(dict1, dict2)
#merged_dict


def merge_dict_probe(data, dict_update, probe_names, probe):
    probe_sub = [item for item in probe_names if probe in item]
    #print(probe_sub)
    dict_tmp = {}
    for item in probe_sub:
        #print(item)
        dict_tmp = merge_dicts(dict_tmp, data[item])
    dict_update[probe] = dict_tmp


def format_percent(x):
    x = str(round(x*100, 2))+"%"
    return x


def run(outdir, sample, json_file, total_amp_umi, total_valid_amp_read):
    with open(json_file, 'r', encoding='utf-8') as file:
        data = json.load(file)
    
    probe_names = sorted(data.keys())

    dict_update = {}
    merge_dict_probe(data, dict_update, probe_names, "16S-V2-R2")
    merge_dict_probe(data, dict_update, probe_names, "16S-V8-R2")
    merge_dict_probe(data, dict_update, probe_names, "515f-R2")
    merge_dict_probe(data, dict_update, probe_names, "P799-R2")
    merge_dict_probe(data, dict_update, probe_names, "V6-F4-R")

    count_amp_dic = dict_update
    count_amp_dic.keys()

    count_amp_list = []
    for amp_name in count_amp_dic.keys():
        UMI_count = 0
        read_count = 0
        #if amp_name in count_amp_dic:
        for cb in count_amp_dic[amp_name]:
            UMI_count += len(count_amp_dic[amp_name][cb])
            read_count += sum(count_amp_dic[amp_name][cb].values())
            #break
        count_amp_list.append({"amp_name": amp_name, "UMI_count": UMI_count, "read_count": read_count})
    df_amp_count = pd.DataFrame(count_amp_list, columns=["amp_name", "read_count", "UMI_count"])

    df_amp_count["read_fraction"] = (df_amp_count["read_count"]/total_valid_amp_read).apply(format_percent)
    df_amp_count["UMI_fraction"] = (df_amp_count["UMI_count"]/total_amp_umi).apply(format_percent)
    df_amp_count.sort_values(by="UMI_count", inplace=True, ascending=False)

    df_amp_count_file = outdir + '/' + sample + '_amp_count_stat.tsv'
    df_amp_count.to_csv(df_amp_count_file, sep="\t", index=False)  



if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', help='outdir', required=True)
    parser.add_argument('--sample', help='sample', required=True)
    parser.add_argument('--json_file', help='json_file', required=True)
    parser.add_argument('--total_amp_umi', help='total_amp_umi', type =int, required=True)
    parser.add_argument('--total_valid_amp_read', help='total_valid_amp_read', type= int, required=True)

    # add version
    parser.add_argument('--version', action='version', version='1.0')
    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        os.system('mkdir -p %s' % (args.outdir))

    run(args.outdir, args.sample, args.json_file, args.total_amp_umi, args.total_valid_amp_read)

