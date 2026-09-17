# PCA 기반 이미지(행렬) 압축 함수
# 입력: x (n x p 행렬), k (근사 rank)
# 출력: 리스트(list) - approx (압축된 행렬), error (Frobenius norm 에러)

pca_compress <- function(x, k) {
  # 1. V = X'X 의 고유분해
  V <- t(x) %*% x
  eig <- eigen(V)
  U <- eig$vectors      # p x p 고유벡터 행렬
  
  # 2. PC 스코어 행렬 Z = XU
  Z <- x %*% U
  
  # 3. rank-k 근사: 처음 k개의 컬럼만 사용
  Uk <- U[, 1:k, drop = FALSE]
  Zk <- Z[, 1:k, drop = FALSE]
  X_k <- Zk %*% t(Uk)
  
  # 4. Frobenius norm 에러
  err <- norm(x - X_k, type = "F")
  
  return(list(approx = X_k, error = err))
}

set.seed(1)
x_test <- matrix(rnorm(20*15), nrow = 20)

result_k5 <- pca_compress(x_test, k = 5)
result_k15 <- pca_compress(x_test, k = 15)

result_k5$error
result_k15$error

library(readr)

# 이미지 CSV → 행렬로 읽기
load_image <- function(path) {
  x <- read_csv(path, show_col_types = FALSE)
  data.matrix(x)
}

# 방향 보정된 이미지 출력 함수
show_image <- function(mat, main = "") {
  image(t(apply(mat, 2, rev)), col = gray.colors(256), axes = FALSE, main = main)
}
img1 <- load_image("image1.csv")

par(mfrow = c(3, 4), mar = c(1,1,2,1))
show_image(img1, main = "Original")

errors_img1 <- numeric(10)
for (k in 1:10) {
  res <- pca_compress(img1, k)
  show_image(res$approx, main = paste0("k = ", k))
  errors_img1[k] <- res$error
}
p <- ncol(img1)
all_errors <- sapply(1:p, function(k) pca_compress(img1, k)$error)

plot(1:p, all_errors, type = "b", pch = 16,
     xlab = "k (rank)", ylab = "Frobenius norm error",
     main = "Approximation Error vs k")
V <- t(img1) %*% img1
eig <- eigen(V)
U <- eig$vectors
Z <- img1 %*% U


plot(U[, 1], type = "l", xlab = "index", ylab = "value",
     main = "First eigenvector (u1)")

plot(Z[, 1], type = "l", xlab = "row index", ylab = "value",
     main = "First PC score (z1)")


# ── 효율적인 에러 곡선 계산 (고유분해를 한 번만 수행) ──
error_curve <- function(x) {
  V <- t(x) %*% x
  eig <- eigen(V)
  U <- eig$vectors
  Z <- x %*% U
  p <- ncol(x)
  
  errors <- numeric(p)
  for (k in 1:p) {
    X_k <- Z[, 1:k, drop = FALSE] %*% t(U[, 1:k, drop = FALSE])
    errors[k] <- norm(x - X_k, type = "F")
  }
  list(errors = errors, U = U, Z = Z)
}

# ── 이미지 하나를 통째로 분석하는 함수 ──
analyze_image <- function(path, image_name = "", k_max = 10) {
  
  img <- load_image(path)
  
  # 1. 원본 + rank 1~10 근사 이미지
  par(mfrow = c(3, 4), mar = c(1,1,2,1))
  show_image(img, main = "Original")
  for (k in 1:k_max) {
    res <- pca_compress(img, k)
    show_image(res$approx, main = paste0("k = ", k))
  }
  
  # 2. 전체 k에 대한 에러 곡선
  par(mfrow = c(1, 1))
  ec <- error_curve(img)
  plot(seq_along(ec$errors), ec$errors, type = "b", pch = 16,
       xlab = "k (rank)", ylab = "Frobenius norm error",
       main = paste0(image_name, ": Approximation Error vs k"))
  
  # 3. 첫 eigenvector / 첫 PC score
  plot(ec$U[, 1], type = "l", xlab = "index", ylab = "value",
       main = paste0(image_name, ": First eigenvector (u1)"))
  
  plot(ec$Z[, 1], type = "l", xlab = "row index", ylab = "value",
       main = paste0(image_name, ": First PC score (z1)"))
  
  invisible(list(img = img, errors = ec$errors, U = ec$U, Z = ec$Z))
}

res1 <- analyze_image("image1.csv", "Image 1")
res2 <- analyze_image("image2.csv", "Image 2")
res3 <- analyze_image("image3.csv", "Image 3")
res4 <- analyze_image("image4.csv", "Image 4")
