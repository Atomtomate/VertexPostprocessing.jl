arr = [1.0 + 1.0im, 2.0 + 2.0im, 3.0+3.0im]

@test VertexPostprocessing.get_conjsymm_f(arr, 0) ≈ 1.0 + 1.0im
@test VertexPostprocessing.get_conjsymm_f(arr, 2) ≈ 3.0 + 3.0im
@test VertexPostprocessing.get_conjsymm_f(arr, -1) ≈ 1.0 - 1.0im
@test VertexPostprocessing.get_conjsymm_f(arr, -1) ≈ 1.0 - 1.0im


