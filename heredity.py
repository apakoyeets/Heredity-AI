import csv
import itertools
import sys

# These are the constant probability values.
PROBS = {
    "gene": {
        2: 0.01,
        1: 0.03,
        0: 0.96
    },
    "trait": {
        2: {True: 0.65, False: 0.35},
        1: {True: 0.56, False: 0.44},
        0: {True: 0.01, False: 0.99}
    },
    "mutation": 0.01
}


def main():
    # Ensure correct usage: exactly one command-line argument (the CSV file).
    if len(sys.argv) != 2:
        sys.exit("Usage: python heredity.py data.csv")
    people = load_data(sys.argv[1])

    # Initialize the probabilities dictionary.
    probabilities = {
        person: {
            "gene": {0: 0, 1: 0, 2: 0},
            "trait": {True: 0, False: 0}
        }
        for person in people
    }

    names = set(people.keys())

    # Loop over all sets of people that might have the trait.
    for have_trait in powerset(names):
        # Filter out cases that conflict with known data.
        fails_evidence = False
        for person in names:
            if people[person]["trait"] is not None:
                if people[person]["trait"] != (person in have_trait):
                    fails_evidence = True
                    break
        if fails_evidence:
            continue

        # Loop over all possible gene assignments.
        for one_gene in powerset(names):
            for two_genes in powerset(names - one_gene):
                # Compute the joint probability for the current assignment.
                p = joint_probability(people, one_gene, two_genes, have_trait)
                # Update probabilities with this joint probability.
                update(probabilities, one_gene, two_genes, have_trait, p)

    # Normalize the results so each probability distribution sums to 1.
    normalize(probabilities)

    # Print results.
    for person in people:
        print(f"{person}:")
        for field in probabilities[person]:
            print(f"  {field.capitalize()}:")
            for value in sorted(probabilities[person][field]):
                p = probabilities[person][field][value]
                print(f"    {value}: {p:.4f}")


def load_data(filename):
    """
    Read a CSV file and return a dictionary where the key is a person's name
    and the value is a dictionary containing:
      name, mother, father, trait
    Empty mother/father fields are converted to None.
    trait is True if "1", False if "0", and None if blank.
    """
    data = {}
    with open(filename) as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row["name"]
            data[name] = {
                "name": name,
                "mother": row["mother"] or None,
                "father": row["father"] or None,
                "trait": (True if row["trait"] == "1" else
                          False if row["trait"] == "0" else None)
            }
    return data


def powerset(s):
    """
    Return a list of all possible subsets of set s.
    Each subset is represented as a Python set.
    """
    s = list(s)
    return [
        set(combo)
        for r in range(len(s) + 1)
        for combo in itertools.combinations(s, r)
    ]


def joint_probability(people, one_gene, two_genes, have_trait):
    """
    Compute and return the joint probability for the given assignment.
    For each person, the probability is:
      - If no parental data:
            Use the unconditional probability in PROBS["gene"]
      - Else, compute the probability based on parent's gene passing:
            For each parent:
              * If parent has 2 copies: probability of passing is (1 - mutation)
              * If parent has 1 copy: 0.5
              * If parent has 0 copies: mutation probability
      - Multiply by the probability of the trait given the gene count.
    """
    probability = 1
    for person in people:
        # Determine the number of genes for this person.
        if person in two_genes:
            genes = 2
        elif person in one_gene:
            genes = 1
        else:
            genes = 0

        # Check if person has the trait.
        has_trait = person in have_trait

        mother = people[person]["mother"]
        father = people[person]["father"]

        if mother is None or father is None:
            # Use unconditional probability.
            p_gene = PROBS["gene"][genes]
        else:
            # Calculate the probability of the gene donation from each parent.
            def pass_prob(parent):
                if parent in two_genes:
                    return 1 - PROBS["mutation"]
                elif parent in one_gene:
                    return 0.5
                else:
                    return PROBS["mutation"]
            mom_prob = pass_prob(mother)
            dad_prob = pass_prob(father)

            if genes == 2:
                p_gene = mom_prob * dad_prob
            elif genes == 1:
                p_gene = mom_prob * (1 - dad_prob) + (1 - mom_prob) * dad_prob
            else:
                p_gene = (1 - mom_prob) * (1 - dad_prob)

        # Multiply in the probability for the trait.
        p_trait = PROBS["trait"][genes][has_trait]
        probability *= p_gene * p_trait

    return probability


def update(probabilities, one_gene, two_genes, have_trait, p):
    """
    Add the joint probability p to the appropriate probabilities for each person.
    """
    for person in probabilities:
        if person in two_genes:
            genes = 2
        elif person in one_gene:
            genes = 1
        else:
            genes = 0
        probabilities[person]["gene"][genes] += p
        probabilities[person]["trait"][person in have_trait] += p


def normalize(probabilities):
    """
    Update the probabilities so that each distribution (gene and trait)
    sums to 1.
    """
    for person in probabilities:
        # Normalize gene probabilities.
        gene_total = sum(probabilities[person]["gene"].values())
        for gene in probabilities[person]["gene"]:
            probabilities[person]["gene"][gene] /= gene_total

        # Normalize trait probabilities.
        trait_total = sum(probabilities[person]["trait"].values())
        for trait in probabilities[person]["trait"]:
            probabilities[person]["trait"][trait] /= trait_total


if __name__ == "__main__":
    main()
