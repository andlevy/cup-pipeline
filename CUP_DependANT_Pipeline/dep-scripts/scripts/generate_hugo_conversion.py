import pandas as pd
import argparse
import json


SYMBOL = "symbol"
PREV = "prev_symbol"
ALIAS = "alias_symbol"


def expand_symbols(df, col_name):
    """Expands a pipe delimited column into multiple rows, e.g. NCRNA00181|A1BGAS|A1BG-AS ->
    NCRNA00181 -
    A1BGAS     -
    A1BG-AS    -"""
    return (
        df[[col_name, SYMBOL]]
        .dropna()
        .assign(**{col_name: lambda x: x[col_name].astype(str).str.split(r"\|")})
        .explode(col_name)
    )


def get_updated_symbols(symbol_list, official_list, prev_mapping_dict, alias_mapping_dict):
    """
    Given a list of gene symbols, return:
    - resolved symbols (official or successfully mapped),
    - unmapped symbols (could not be resolved),
    - original list with symbols updated where possible (same order as input)
    """
    symbol_map = {}
    unmapped = []

    for symbol in symbol_list:
        if symbol in official_list:
            true = symbol
        elif symbol in prev_mapping_dict:
            true = prev_mapping_dict[symbol]
        elif symbol in alias_mapping_dict:
            true = alias_mapping_dict[symbol]
        else:
            true = False
            unmapped.append(symbol)  # Keep original if it couldn't be updated

        if true:
            symbol_map[true] = true
            prev = {x for x in prev_mapping_dict if prev_mapping_dict[x] == true and x not in official_list}
            alias = {x for x in alias_mapping_dict if alias_mapping_dict[x] == true and x not in official_list}
            for x in prev.union(alias):
                symbol_map[x] = true

    return symbol_map, unmapped


def map_from_old_to_symbol(hgnc_df, old_column=PREV):
    old_df = hgnc_df[hgnc_df[old_column].notna()]
    expanded_df = expand_symbols(old_df, old_column)
    old_to_official = (
        expanded_df
        .set_index(old_column)[SYMBOL]
        .to_dict()
    )
    return old_to_official


def generate_hugo_conversion(original_hgnc, depmap_mutations, output):
    hgnc = pd.read_csv(original_hgnc, sep="\t", low_memory=False)
    mutation_df = pd.read_csv(depmap_mutations, low_memory=False)

    alias_symbol_to_official = map_from_old_to_symbol(hgnc, ALIAS)
    prev_symbol_to_official = map_from_old_to_symbol(hgnc, PREV)

    official_symbols = set(hgnc[hgnc[SYMBOL].notna()][SYMBOL].tolist())

    mutation_gene_symbols = mutation_df['HugoSymbol'].unique().tolist()

    symbol_mapping, unmapped = get_updated_symbols(mutation_gene_symbols,
                                                   official_symbols,
                                                   prev_symbol_to_official,
                                                   alias_symbol_to_official)

    print(f"{len(unmapped)} depmap genes did not have appropriate hugo symbols")
    with open(output, 'w') as f:
        json.dump(symbol_mapping, f)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert hgnc file into a hugo conversion file for depmap genes")
    parser.add_argument('--original', type=str, help='hgnc file')
    parser.add_argument('--depmap', type=str, help='Omics Somatic Mutation file')
    parser.add_argument('--output', type=str, help='Where to write conversion dictionary to')

    args = parser.parse_args()
    generate_hugo_conversion(args.original, args.depmap, args.output)
