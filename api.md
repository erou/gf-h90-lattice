# API Flint

    solve_h90
        INPUT: ζ a n-th root of unity, 
        OUTPUT: x a solution of H90 for ζ, with x^n = the right constant
    
    norm_r
        INPUT: x an element of the Kummer algebra, (a, b) two integers
        OUTPUT: the right norm N_{b/a} of x
    
    first_coordinate
        INPUT: x an element of the Kummer algebra, ζ a n-th root of unity
        OUTPUT: the first coordinate of x in the base 1⊗ζ
    
    matrices (fq_embed_matrices)
        INPUT: x, y elements of two finite fields
        OUTPUT: M, N two matrices corresponding to the map x|->y and its section

## Julia

In Julia there will be global variables
 - `EMBEDDINGS` containing already computed embeddings
 - `H90_ELEMENTS` containing already computed solutions of H90
 - `ROOTS_UNITY` containing the chosen roots of unity

A function `embed` that will take care of everything:

    embed
        INPUT: k, K two finite fields
        OUTPUT: an embedding from k to K
    1. construct the Kummer algebras A and B associated to k and K using the roots in
       ROOTS_UNITY
    2. Check if there is an embedding from k to K in EMBEDDINGS, if so, return
       it
    3. Check if there is a solution alpha of H90 in A, if not compute it using
       solve_h90 and save it in H90_ELEMENTS
    4. Check if there is a solution beta of H90 in B, if not compute it using
       solve_h90 and save it in H90_ELEMENTS
    5. Compute f(beta), where f is a function using an exponantiation, or a norm, depending
       on the degrees of k and K
    6. Compute the first coordinates u and v of alpha and f(beta)
    7. Compute the embedding u|->v thanks to the function matrice
