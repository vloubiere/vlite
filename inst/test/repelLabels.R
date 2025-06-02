sel <- readRDS("/groups/stark/vloubiere/projects/epigenetic_cancer/db/model/ATAC_FC_PH29_PHD11_motif_lm.rds")
x <- sel$`t value.d11`
y <- sel$`t value.29`
labels <- sel$name
cex <- 1
plot(x, y)

dat <- data.table(x,
                  y,
                  labels,
                  cex)
dat[, width:= strwidth(labels, cex= cex)]
dat[, height:= strheight(labels, cex= cex)]
dat[, x1:= x-width/2]
dat[, y1:= y-height/2]
dat[, x2:= x+width/2]
dat[, y2:= y+height/2]
# rect(dat$x1, dat$y1, dat$x2, dat$y2)
# text(x, y, labels)

# Points and boxes DF ----
posdf <- data.frame(x = x, y = y)
boxdf <- as.data.frame(dat[, x1:y2])

# Compute positions boxes
moved <- latticetools::repel_boxes(data_points = as.matrix(posdf),
                                   point_size = rep(0.1, nrow(posdf)),# npc unit
                                   point_padding_x = 0.01,
                                   point_padding_y = 0.01,
                                   boxes = as.matrix(boxdf),
                                   xlim = par("usr")[c(1,2)],
                                   ylim = par("usr")[c(3,4)],
                                   hjust = rep(0, nrow(posdf)),
                                   vjust = rep(0, nrow(posdf)),
                                   force_push = 1e-3,
                                   force_pull = 1e-3,
                                   max_time = 3,
                                   max_overlaps = 10,
                                   max_iter = 100000,
                                   direction = "both",
                                   verbose = FALSE)
setnames(moved, c("x", "y"), c("x_moved", "y_moved"))
dat <- cbind(dat, moved)
dat[, x1:= x1+(x_moved-x)]
dat[, x2:= x2+(x_moved-x)]
dat[, y1:= y1+(y_moved-y)]
dat[, y2:= y2+(y_moved-y)]
# rect(dat$x1, dat$y1, dat$x2, dat$y2)
# text(dat$moved_x, dat$moved_y, labels)

# Compute clipping masks for segments ----
dat[, clip_x1:= par("usr")[1]]
dat[, clip_x2:= par("usr")[2]]
dat[, clip_y1:= par("usr")[3]]
dat[, clip_y2:= par("usr")[4]]
dat[, max.y.diff:= ((y2-y1)/2)*(abs(x-x_moved)/((x2-x1)/2))]
dat[(y-y_moved)>=max.y.diff, clip_y1:= y2]
dat[(y-y_moved)<=(-max.y.diff), clip_y2:= y1]
dat[abs(y-y_moved)<=max.y.diff & x>=x2, clip_x1:= x2]
dat[abs(y-y_moved)<=max.y.diff & x<=x1, clip_x2:= x1]

# Plot segments ----
dat[, seg_length.x:= min(abs(c(x1, x2)-x)), .(x1, x2, x)] # Check if segments are long enough
dat[, seg_length.y:= min(abs(c(y1, y2)-y)), .(y1, y2, y)]
dat[, seg_length:= seg_length.x>strwidth("M", cex= label.cex/3) | seg_length.y>strheight("M", cex= label.cex/3)]
if(any(dat$seg_length))
{
  dat[(seg_length),
      {
        # clip(clip_x1[1], clip_x2[1], clip_y1[1], clip_y2[1])
        segments(x_moved[1], y_moved[1], x[1], y[1], col= seg.col[1])
      }, .(x_moved, y_moved, x, y, clip_x1, clip_x2, clip_y1, clip_y2, seg.col)]
}
# clip(par("usr")[1], par("usr")[2], par("usr")[3], par("usr")[4])
# Plot rectangles
if(rect.draw)
  with(dat, rect(x1, y1, x2, y2, col = rect.col, border= NA))
# Plot labels
with(dat, text(x, y, labels = label, col = label.col, cex = label.cex))
# Return object
invisible(dat)

# return coordinates
return(dat)
