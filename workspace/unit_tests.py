from unittest import TestCase
from subprocess import call
from create_database import pattern_filter, taxonomy_information


class TestDatabase(TestCase):
    "Tests on methods to index genomes for database"

    def test_apply_pattern(self):
        "Tests if a pattern is applied"
        self.assertEqual(pattern_filter("ATCAG", [1, 1, 0, 1, 1]), "ATAG")

    def test_extract_taxo(self):
        "Tests if a taxonomy is correctly extracted"
        self.assertEqual(
            taxonomy_information(
                "path/to/some/file/Bacteria_RealBacteria_TrueBacteria_SomeBacteria_AmazingBacteria.fna"
            ),
            {
                'domain': 'Bacteria',
                'phylum': 'RealBacteria',
                'group': 'TrueBacteria',
                'order': 'SomeBacteria',
                'family': 'AmazingBacteria'
            }
        )


if __name__ == "__main__":
    call("python -m unittest -v unit_tests.py", shell=True)
