
ww <- seq(1, 5)
a <- c(0.5, 0.55, 0.55, 0.6, 0.6)

png(filename = "icon.png",
    width = 250, height = 100, units = "px", pointsize = 20,
    bg = "white", res = NA, family = "", restoreConsole = TRUE,
    type = c("windows", "cairo", "cairo-png"))
par(mar = c(0, 0, 0, 0), bg = NA)
plot(ww, a, type = "o", lwd = 3, ylim = c(0.4, 0.8), xlim = c(1, 7), bty = "n", xaxt = "n", yaxt = "n", ylab = "", xlab = "", col = "white")

abline(h = 0.75, lty = 2, lwd = 2, col = "#F8CD00")
abline(h = 0.45, lty = 2, lwd = 2, col = "#F8CD00")

arrows(x0 = 5, y0 = 0.6, x1 = 6, y1 = 0.7, lwd = 3, col = "white", length = 0.1)

dev.off()

library(hexSticker)

s <- sticker("icon.png",
             package="sequential.pops", p_size=15, h_fill = "#2C5698", h_color = "#F8CD00",s_x=1, s_y=.8, s_width=1, s_height=.4,
             filename="sticker.png")

