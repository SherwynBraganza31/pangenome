import os, gffutils, pandas as pd, json

def updateDatabase(db: pd.DataFrame, attributes: dict, organism: str) -> bool:
    """
      Takes in an attributes dictionary and adds data from it as a row
      in the db DataFrame.
      If the sub-attributes that its looking for doesn't exist,
      it skips the update and returns a False; else True.

      Attribute Checks:
      - Checks to see the dict contains data for the db cols (col_keys)
      - checks to see that inference is a list of 2 items -- The UniProtKB item
      is only present if inferece as the second element:

      All elements of the dictionary are lists of strings.

      @params
          db: pandas.DataFrame
              DataFrame that represents the main 'database' of cherry
              picked data
          attributes: dict
              Dictionary corresponding to the attributes of a row of the
              .gff file
          organims: str
              The organism from which the gene came from

        @returns
          True:  if attributes dictionary passes all the checks and updates
          False: if it fails the checks

    """
    col_keys = ['ID', 'Name', 'inference', 'product']

    if 'ID' in attributes.keys():
        ID, Name, UniProtKB, product = None, None, None, None

        # split the string by :, the last element contains the UniProtKB
        for substring in attributes['inference']:
            if 'UniProtKB' in substring:
                UniProtKB = substring.split(':')[-1]
                break

        try:
            ID = attributes['ID'][0]
        except KeyError:
            ID = attributes['locus-tag']

        try:
            Name = attributes['Name'][0]
        except KeyError:
            pass

        try:
            product = attributes['product'][0]
        except KeyError:
            pass

        db.loc[len(db.index)] = [ID, Name, UniProtKB, product, organism]
        return True
    else:
        return False


if __name__ == "__main__":
    parent_path = input("Enter the parent directory")
    parent_path = parent_path if "/" in parent_path[-1] else parent_path + "/"
    gff_path = parent_path + "ppan_dataset/GFF/"
    gff_file_list = os.listdir(gff_path)
    os.chdir(gff_path)

    gene_count = {}

    out_df = pd.DataFrame(columns=['Gene ID', 'Gene Symbol', 'UniProtKB',
                                   'Annotation', 'Organism'])

    for file in gff_file_list:
        count = 0
        if '.gff' not in file:
            continue
        # create a sqlite3 db of the file in memory
        curr_file = gffutils.create_db(file, dbfn=':memory:',
                                       force=True, keep_order=True, merge_strategy='merge',
                                       sort_attribute_values=True)
        # sql stmt to select attributes of all rows
        for row in curr_file.execute('SELECT attributes FROM features'):
            # parse string attributes into dict
            attribs = json.loads(row['attributes'])
            res = updateDatabase(out_df, attribs, file.split('.')[0])
            count+=1
        gene_count.update({file:count})

    print(gene_count)
    out_df.set_index('Gene ID', inplace=True)
    out_df.to_csv(parent_path + "trimmed_matrix_files/gff_sequencing.csv")