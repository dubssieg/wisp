from .data_manipulation import species_list
from .data_manipulation import species_map
from .data_manipulation import load_xgboost_data
from .data_manipulation import write_xgboost_data
from .data_manipulation import check_if_database_exists
from .data_manipulation import check_if_merged_database_exists
from .data_manipulation import check_if_merged_model_exists
from .data_manipulation import check_if_model_exists
from .data_manipulation import load_mapping

from .kmers_coders import encode_kmer_4
from .kmers_coders import decode_kmer_4
from .kmers_coders import recode_kmer_4
from .kmers_coders import kmer_indexing_brut
from .kmers_coders import read_and_its_compl
from .kmers_coders import optimal_splitting
from .kmers_coders import reverse_comp
from .kmers_coders import splitting_generator
from .kmers_coders import counter_ultrafast
from .kmers_coders import encoder

from .parameters_init import load_json
