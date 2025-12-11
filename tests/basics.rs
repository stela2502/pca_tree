use pca_tree::PcaTree;

#[test]
fn small_sequences_produce_tree() {
    let seqs = vec![
        "AAAA".to_string(),
        "AAAT".to_string(),
        "AATT".to_string(),
        "TTTT".to_string(),
    ];

    let mut model = PcaTree::new(seqs, 2);
    model.fit().unwrap();

    let coords = model.coords().unwrap();

    // Should have N rows
    assert_eq!(coords.nrows(), 4);

    // Should have k columns
    assert_eq!(coords.ncols(), 2);

    // Tree must have N-1 edges
    assert_eq!(model.tree().len(), 3);

    // All nodes should appear at least once
    let mut seen = vec![false; 4];
    for (p, c, _) in model.tree() {
        seen[*p] = true;
        seen[*c] = true;
    }

    assert!(seen.iter().all(|x| *x));
}

