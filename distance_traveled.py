import utils


suffixes = {"genotypes": "_genotypes.csv", "fish_used": "_fish.txt"}
all_files = utils.find_data("distance_traveled", suffixes)

main_files = all_files["main_files"]
wildtype_files = all_files["wildtype_files"]

def analyze(files):
    for prefix_name in files.keys():
        group = files[prefix_name]
        genotypes = utils.load_genotypes(group["genotypes"], group["fish_used"], "down")
        print(genotypes)

analyze(main_files)
