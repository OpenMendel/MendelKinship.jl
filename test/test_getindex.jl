using BenchmarkTools

function correspond(x::Vector{Int}, y::Vector{Int})
  (m, n) = (length(x), length(y))
  xperm = sortperm(x)
  yperm = sortperm(y)
  x_to_y = zeros(Int, m)
  # y_to_x = zeros(Int, n)
  (i, j) = (1, 1)
  done = false
  while !done
    ii = xperm[i]
    jj = yperm[j]
    if x[ii] == y[jj]
      x_to_y[ii] = jj
      # y_to_x[jj] = ii
      (i, j) = (i + 1, j + 1)
    elseif x[ii] < y[jj]
      i = i + 1
    elseif y[jj] < x[ii]
      j = j + 1
    end
    done = i > m #|| j > n
  end
  # return (x_to_y, y_to_x)
  return x_to_y
end

x = collect(1:10^5)
y = collect(1:10^5)
shuffle!(x)
shuffle!(y)

@benchmark correspond(x, y)
@benchmark indexin(x, y)
[indexin(x, y), indexin(y, x)]





