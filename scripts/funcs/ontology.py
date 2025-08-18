from owlready2 import *
import pandas as pd

def load_ontology(ontology, namespace):
    """
    Load ontology
    """

    onto = get_ontology(ontology)
    namespace = get_namespace(namespace)

    if onto.loaded:
        print("Ontology is already loaded.")
    else:
        onto.load()
        print("Ontology was not loaded. Now loaded.")

    return onto, namespace

def get_term_from_label(namespace, label, ontology, prefix = 'obo.'):
    """
    Get ontology term from label
    """
    result = namespace.search(label=label)
    if result:
        id = result[0]  # Returns the IRI of the first matching result
        id = str(id).replace(prefix,'')  # Remove 'obo.' from the IRI
        id = id.replace(f'{ontology}.', '')  # Remove the ontology IRI from the IRI
        return id
    else:
        return None  # Return None if no result found
        
def search_term(onto, namespace,search_string):
    # Case-insensitive search through all classes
    matching_classes = [
        cls for cls in onto.classes()
        if any(search_string.lower() in label.lower() for label in cls.label)
    ]

    # Print results
    for cls in matching_classes:
        print(f"IRI: {cls.iri}")
        print(f"Label(s): {cls.label}")
        print(f"ID: {cls.name}")

    #create a DataFrame to store results
    data = {
        'GO_Label': [cls.label[0] for cls in matching_classes],
        'GO_ID': [cls.name.replace('_', ':') for cls in matching_classes]
    }
    df = pd.DataFrame(data)
    print(df)
    return df
