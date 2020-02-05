from typing import List, Dict, Union
import gzip
import json


def print_header_row():
    print(r'\begin{tabular}{lll}')
    print(r'\vline')
    columns: List[str] = ['Data set', 'Prior', 'Top Covariates (Covariate:Inclusion probability)']
    print('&'.join(columns) + r'\\')
    print(r'\vline')


def print_table_rows(top_inclusions: List[Dict[str, Union[str, Dict[str, List[Union[float, int]]]]]]):
    for idx, row in enumerate(top_inclusions):
        data_set: str = str(row['data_sets'])
        prior: str = str(row['priors'])
        inclusion_probs: List[float] = row['top_inclusions']['vinclusion_probs']
        top_columns: List[int] = row['top_inclusions']['top_columns']
        print_table_row(data_set, prior, inclusion_probs, top_columns)
        if idx != len(top_inclusions) - 1:
            print(r'\\')
        else:
            print('')


def print_table_row(data_set: str, prior: str, inclusion_probs: List[float], top_columns: List[int]):
    # columns_inclusions: str = ''
    assert len(inclusion_probs) == len(top_columns)
    # FIXME - This isn't really good enough. You should make a table within a
    # table. Like this https://tex.stackexchange.com/questions/7958/how-to-nest-tables
    data_set_sanitised: str = data_set.replace('_', ' ') 
    print(f'{data_set_sanitised} & {prior} & ', end='')
    column_format: str = 'l' * len(top_columns)
    print(r'\begin{tabular}{%s}' % column_format)
    for idx in range(len(top_columns)):
      print(top_columns[idx], end='')
      if idx < len(top_columns) - 1: print('&', end='')
    print(r'\\')
    for idx in range(len(inclusion_probs)):
      print(inclusion_probs[idx], end='')
      if idx < len(inclusion_probs) - 1: print('&', end='')
    print(r'\\')
    print(r'\end{tabular}')
    # for idx in range(len(top_columns)):
    #     columns_inclusions += f'{top_columns[idx]}:{inclusion_probs[idx]} '
    # print(data_set, prior, columns_inclusions, sep='&', end='')


def print_footer_row():
    print(r'\vline')
    print(r'\end{tabular}')


def main():
    top_inclusions: List[Dict[str, Union[str, Dict[str, List[Union[float, int]]]]]]
    top_inclusions = json.load(gzip.open('top_inclusions_df.json.gz'))
    print_header_row()
    print_table_rows(top_inclusions)
    print_footer_row()

if __name__ == "__main__":
    main()
