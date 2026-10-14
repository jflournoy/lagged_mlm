# Path diagram of the data-generating model, drawn from the generator's own
# lavaan syntax (sim_factorial_lib.R's make_generator_text()) so the picture
# can't drift from what is simulated. Paths the factorial varies are orange,
# the y -> x path (always 0) is dashed, and the trait correlation is drawn as
# a bracket outside the waves.

dgp_diagram <- function(W = 12) {
  pt <- lavaan::lavaanify(make_generator_text(list(waves = W, sx = .7, sy = .7, clxy = 0, clyx = 0,
                                                    co = .5, v = 1, tcor = 0)),
                          int.ov.free = FALSE, auto.var = FALSE, auto.fix.first = FALSE)
  pt <- pt[pt$op %in% c("=~", "~", "~~") & !(pt$op == "~~" & pt$lhs == pt$rhs), ]
  e <- data.frame(from = ifelse(pt$op == "~", pt$rhs, pt$lhs), to = ifelse(pt$op == "~", pt$lhs, pt$rhs),
                  arrow = ifelse(pt$op == "~~", "both", "last"), stringsAsFactors = FALSE)
  # Covariances among observed variables are fixed at 0 in the generator.
  e <- e[!(grepl("^[xy][0-9]", e$from) & grepl("^[xy][0-9]", e$to) & e$arrow == "both"), ]
  role <- with(e, ifelse(from %in% c("kappa", "omega") & arrow == "last", "trait loading",
            ifelse(grepl("^[pq]", from) & grepl("^[xy]", to), "within part",
            ifelse(grepl("^p", from) & grepl("^p", to), "x stability",
            ifelse(grepl("^q", from) & grepl("^q", to), "y stability",
            ifelse(grepl("^p", from) & grepl("^q", to) & arrow == "last", "x to y",
            ifelse(grepl("^q", from) & grepl("^p", to) & arrow == "last", "y to x",
            ifelse(arrow == "both" & grepl("^[pq]", from), "within-wave r",
            ifelse(arrow == "both" & from %in% c("kappa", "omega"), "trait r", NA)))))))))
  if (anyNA(role)) stop("unclassified edge in the generator: ", paste(e$from[is.na(role)], e$to[is.na(role)]))
  keep <- role != "trait r"
  e <- e[keep, ]; role <- role[keep]
  e$colour <- ifelse(role %in% c("y stability", "x to y", "within-wave r"), "#d4572a",
              ifelse(role == "y to x", "#9a9a9a", "#555555"))
  e$linetype <- ifelse(role == "y to x", 2, 1)
  e$curvature <- ifelse(role == "within-wave r", 60, NA)
  e$label <- ""
  lay <- matrix(NA_character_, nrow = 6, ncol = W)
  lay[1, W / 2] <- "kappa"; lay[2, ] <- paste0("x", 1:W); lay[3, ] <- paste0("p", 1:W)
  lay[4, ] <- paste0("q", 1:W); lay[5, ] <- paste0("y", 1:W); lay[6, W / 2] <- "omega"
  nodes <- data.frame(name = c("kappa", "omega", paste0("x", 1:W), paste0("y", 1:W), paste0("p", 1:W), paste0("q", 1:W)),
                      shape = c("oval", "oval", rep("rect", 2 * W), rep("oval", 2 * W)), stringsAsFactors = FALSE)
  nodes$label <- c("x trait\n(var 0, .5, 1)", "y trait\n(var 0, .5, 1)", paste0("x", 1:W), paste0("y", 1:W),
                   paste0("x", 1:W, "*"), paste0("y", 1:W, "*"))
  g <- tidySEM::prepare_graph(edges = e, nodes = nodes, layout = lay, rect_width = .9, rect_height = .6,
                              ellipses_width = .9, ellipses_height = .6, spacing_x = 1.4, spacing_y = 1.4,
                              text_size = 3.2, label_size = 3)
  nd <- tidySEM::nodes(g)
  nd$colour <- ifelse(nd$name %in% c("kappa", "omega"), "#d4572a", "#444444")
  tidySEM::nodes(g) <- nd
  k <- nd[nd$name == "kappa", ]; o <- nd[nd$name == "omega", ]
  xr <- max(nd$node_xmax) + .6
  plot(g) +
    ggplot2::annotate("path", x = c(k$node_xmax, xr, xr, o$node_xmax), y = c(k$y, k$y, o$y, o$y),
                      colour = "#d4572a", linewidth = .5,
                      arrow = grid::arrow(ends = "both", length = grid::unit(.12, "inches"), type = "closed")) +
    ggplot2::annotate("text", x = xr + .15, y = (k$y + o$y) / 2, label = "trait r",
                      hjust = 0, size = 3.2, colour = "#d4572a") +
    dgp_key(x0 = xr + 1.4, ytop = k$y, ybottom = o$y)
}

# Key: each path type, its value in the main simulations, and (orange) its
# values in the factorial.
dgp_key <- function(x0, ytop, ybottom) {
  key <- data.frame(
    label = c("x stability: .7", "y stability: .7  (factorial: 0, .3, .7)", "x \u2192 y: 0  (factorial: 0, .1)",
              "y \u2192 x: 0", "within-wave r: .5  (factorial: 0, .5)", "trait r: 0  (factorial: 0, .3)",
              "trait variance: 0, .5, 1"),
    colour = c("#555555", "#d4572a", "#d4572a", "#9a9a9a", "#d4572a", "#d4572a", "#d4572a"),
    linetype = c(1, 1, 1, 2, 1, 1, 0), stringsAsFactors = FALSE)
  y <- seq(ytop - .4, ybottom + .4, length.out = nrow(key))
  list(ggplot2::annotate("text", x = x0, y = ytop + .35, label = "Paths (main simulations; factorial)",
                         hjust = 0, size = 3.4, fontface = "bold", colour = "#333333"),
       ggplot2::annotate("segment", x = x0, xend = x0 + .7, y = y[key$linetype > 0], yend = y[key$linetype > 0],
                         colour = key$colour[key$linetype > 0], linetype = key$linetype[key$linetype > 0], linewidth = .6,
                         arrow = grid::arrow(length = grid::unit(.08, "inches"), type = "closed")),
       ggplot2::annotate("text", x = x0 + .85, y = y, label = key$label, hjust = 0, size = 3.2, colour = key$colour),
       ggplot2::expand_limits(x = x0 + 6.2))
}
