# datasets=("syn-uni/default" "syn-gau/default" "syn-zipf/default" "yeast" "hprd" "human" "wordnet" "dblp" "youtube" "patents")

# datasets=("syn-uni/default" "syn-gau/default" "syn-zipf/default" "yeast" "hprd" "dblp" "youtube" "patents")
# datasets=("syn-uni/default" "syn-gau/default" "syn-zipf/default")
datasets=("yeast" "hprd" "dblp" "youtube" "patents")
for dataset in "${datasets[@]}"; do
    d_arg="../Datasets/${dataset}"

    ./build/NSM -d $d_arg -q cyclic/8-3
done

# echo "Parameter Tuning \n hop num"

# datasets=("syn-uni/default" "syn-gau/default" "syn-zipf/default")
# for dataset in "${datasets[@]}"; do
#     d_arg="../Datasets/${dataset}"

#     ./build/NSM -d $d_arg -q cyclic/8-3 -k 1
#     ./build/NSM -d $d_arg -q cyclic/8-3 -k 3
#     ./build/NSM -d $d_arg -q cyclic/8-3 -k 4
# done

# echo "ratio"
# for dataset in "${datasets[@]}"; do
#     d_arg="../Datasets/${dataset}"

#     ./build/NSM -d $d_arg -q cyclic/8-3 -a 10 -b 1
#     ./build/NSM -d $d_arg -q cyclic/8-3 -a 100 -b 0.1
#     ./build/NSM -d $d_arg -q cyclic/8-3 -a 100 -b 0.01
#     ./build/NSM -d $d_arg -q cyclic/8-3 -a 10000 -b 0.01
# done

# datasets=("syn-uni/default" "syn-gau/default" "syn-zipf/default")
# for dataset in "${datasets[@]}"; do
#     d_arg="../Datasets/${dataset}"

#     ./build/NSM -d $d_arg -q cyclic/5-3
#     ./build/NSM -d $d_arg -q cyclic/6-3
#     ./build/NSM -d $d_arg -q cyclic/10-3
#     ./build/NSM -d $d_arg -q cyclic/12-3
# done

# echo ""
# echo ""
# echo ""

# datasets=("syn-uni/default" "syn-gau/default" "syn-zipf/default")
# for dataset in "${datasets[@]}"; do
#     d_arg="../Datasets/${dataset}"

#     ./build/NSM -d $d_arg -q cyclic/8-2
#     ./build/NSM -d $d_arg -q cyclic/8-4
# done

# distributions=("syn-uni")
# datasets=("size-10K" "size-30K" "default" "size-80K" "size-100K" "size-500K" "size-1M" "size-5M" "size-10M")
# # distributions=("syn-uni" "syn-gau" "syn-zipf")
# # datasets=("deg-3" "deg-4" "deg-6" "deg-7" "size-10K" "size-30K" "size-80K" "size-100K" "size-500K" "size-1M" "size-5M" "size-10M")
# # datasets=("label-5" "label-10" "label-20" "label-25")

# for distribution in "${distributions[@]}"; do
#     for dataset in "${datasets[@]}"; do
#         d_arg="../Datasets/${distribution}/${dataset}"

#         ./build/NSM -d "$d_arg" -q "cyclic/8-3"

#         echo ""
#     done

#     echo -e "\n\n\n"
# done