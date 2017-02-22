lapply(seq_along(b$csize)[b$csize > 1], function(x) 
  V(a)$name[b$membership %in% x])

table(names(V(ESC)) %in% names(V(PRK)))
