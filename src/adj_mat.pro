FUNCTION adj_mat, matrix
  return, [[matrix[1, 1],-matrix[1, 0]],[-matrix[0, 1], matrix[0, 0]]]
END