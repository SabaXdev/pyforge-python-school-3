import pytest
from main import substructure_search


@pytest.mark.parametrize("smiles, substructure, expected", [
    # Empty List
    ([], "CCO", []),
    # Smile with no matches
    (["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "CH4", []),
    # Basic test cases
    (["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"],
     "c1ccccc1", ["c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"]),
    (["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"],
     "CO", ['CCO', 'CC(=O)O', 'CC(=O)Oc1ccccc1C(=O)O']),
    # Edge case: very simple SMILES strings
    (["C", "O", "CC", "CO"], "C", ["C", "CC", "CO"]),
])
def test_substructure_search(smiles, substructure, expected):
    assert substructure_search(smiles, substructure) == expected
