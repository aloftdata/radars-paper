Case Study weather radar
================

## Setup

``` r
require(bioRad)
require(getRad)
require(lubridate)
require(dplyr)
require(sf)
require(ggplot2)
require(stars)
```

In this case study we zoom in to the morning of the 19th of October
2022. We select both radar in the Netherlands. For a quick scan of the
data the crow visualizer
[aloft](https://crow.aloftdata.eu#/?radar=bejab&date=2022-10-19&interval=24&timedisplay=radarLocal&vpiMode=mtr&vpColorScheme=viridis&lang=en)
repository can be consulted.

``` r
time <- as_datetime("2022-10-19 8:00:00")
radarCode <- c("nldhl", "nlhrw")
radarLocation <- get_weather_radars() |>
  dplyr::filter(radar %in% radarCode) |>
  st_transform(crs = 28992)
radarGrid <- raster::raster(st_buffer(radarLocation, units::set_units(105, "km")),
  res = 500
)
```

On this night high densities of numbers were counted along the dutch
coast, as seen in the following video on youtube.

[![Watch the
video](https://img.youtube.com/vi/8x4Dos7GbvI/default.jpg)](https://youtu.be/8x4Dos7GbvI)

For visualizations we load a polygon map, in this case from the
`cbsodataR` package. This is a slightly nicer map compared to
`rnaturalearth` which globally provides reasonable maps.

``` r
require(cbsodataR)
```

    ## Loading required package: cbsodataR

``` r
map <- cbs_get_sf("provincie", year = 2022)
```

# Mornings activity

To explore the mornings activbity we can use the vpts.

``` r
vpts <- get_vpts(radar = radarCode, as.Date(time))
plot(vpts$nlhrw, main="Herwijnen")
```

![](case_study_3_reproducible_files/figure-gfm/vpts-1.png)<!-- -->

``` r
plot(vpts$nldhl, main="Den Helder")
```

![](case_study_3_reproducible_files/figure-gfm/vpts-2.png)<!-- -->

# Static map

The first step is to download radar data and calculate vertical
distributions from the data

``` r
# Set names is to ensure the resulting list is named according to the radar code
pvol_list <- setNames(
  get_pvol(
    radarCode, time,
    elev_max = 89
  ), radarCode
)
vp_list <- purrr::map(pvol_list, calculate_vp)
```

    ## Running vol2birdSetUp

    ## Warning: radial velocities will be dealiased...

    ## Running vol2birdSetUp

    ## Warning: radial velocities will be dealiased...

The next step is to create for each radar separately a vertically
integrated estimate of density. Here we convert all cases where strong
corrections have been applied to `NA`

``` r
ppi_strs <- purrr::map2(pvol_list, vp_list, integrate_to_ppi, raster = radarGrid) |>
  purrr::map(purrr::pluck, "data") |>
  purrr::map(stars::st_as_stars) |>
  c(list(along = "radar")) |>
  do.call(what = "c") |>
  mutate(VID = dplyr::if_else(R > 10, NA, VID))
```

    ## Loading required namespace: raster

To combine both radars we create a dataset where we identify which radar
on which location as been minimally corrected.

``` r
minR <- st_apply(select(ppi_strs, R), 1:2, function(x) {
  if (all(is.na(x))) {
    x
  } else {
    x == min(x, na.rm = T)
  }
}, keep = T) |>
  setNames("minR") |>
  aperm(c(2, 3, 1))
```

By multiplying the VID with the minimal R value (0/1)and summing the
values for each location we can get the best VID estimate for each
location

``` r
ppi_minR <- c(ppi_strs, minR) |>
  mutate(VID = VID * minR) |>
  st_apply(c("x", "y"), sum, na.rm = T)
```

``` r
ggplot() +
  geom_stars(data = ppi_minR, aes(fill = VID), color = NA, sf = T) +
  scale_fill_viridis_c(
    name = "Bird density [1/km^2]",
    option = "turbo",
    limit = c(2, 800),
    transform = "log10"
  ) +
  geom_sf(data = map, fill = NA, color = "black", size = .3) +
  geom_sf(data=radarLocation)+
  ggspatial::annotation_scale() +
  coord_sf(
    expand = F,
    xlim = st_bbox(radarGrid)[c("xmin", "xmax")],
    ylim = st_bbox(radarGrid)[c("ymin", "ymax")]
  )
```

    ## Loading required namespace: raster

    ## Warning in scale_fill_viridis_c(name = "Bird density [1/km^2]", option =
    ## "turbo", : log-10 transformation introduced infinite values.

![](case_study_3_reproducible_files/figure-gfm/exampleRplot-1.png)<!-- -->

In this map we see high density of birds along the coast this
corresponds to [citizen science
observations](https://trektellen.nl/maps/species/-2/0/0/20221019/20221019/1/0)
along the coast. This is mainly due to an exceptional number of redwing
passing by as noted by
[SOVON](https://www.vogelbescherming.nl/actueel/bericht/koperwiekentrek).

# Animation

For the animation we first generate a vector with all timestamps we want
to visualize

``` r
animation_interval <- interval(
  time - lubridate::minutes(180),
  time + lubridate::minutes(20)
)
animation_timestamps <- seq(
  int_start(animation_interval),
  int_end(animation_interval),
  "5 mins"
)
```

The next step is a function to translate each timestamps in a dataset
for visualization. This function is based on the example before, however
extra attention has been given to make the code memory efficient. For
larger timeseries you might want to implement caching of intermediate
results.

``` r
vid_for_timestamp <- function(x) {
  pvol_list <- getRad::get_pvol(
    radarCode,
    x,
    param = c("DBZH", "VRADH", "RHOHV"),
    elev_max = 89
  ) |>
    setNames(radarCode)
  vp_list <- purrr::map(pvol_list, bioRad::calculate_vp, warnings = F)
  ppi_list <- purrr::map2(
    pvol_list,
    vp_list,
    bioRad::integrate_to_ppi,
    raster = radarGrid
  ) |>
    purrr::map(purrr::pluck, "data") |>
    purrr::map(stars::st_as_stars) |>
    c(list(along = "radar")) |>
    do.call(what = "c") |>
    mutate(VID = dplyr::if_else(R > 10, NA, VID)) |>
    select(R, VID)
  minR <- st_apply(
    select(ppi_list, R),
    1:2,
    function(x) {
      if (all(is.na(x))) {
        x
      } else {
        x == min(x, na.rm = T)
      }
    },
    keep = T
  ) |>
    setNames("minR") |>
    aperm(c(2, 3, 1))
  ppi_minR <- c(ppi_list, minR) |>
    mutate(VID = VID * minR) |>
    st_apply(c("x", "y"), sum, na.rm = T) |>
    select(VID) |>
    sf::st_as_sf(as_points = T) |>
    tibble::add_column(time = x) |>
    filter(!is.na(VID)) |>
    mutate(
      X = sf::st_coordinates(geometry)[, "X"],
      Y = sf::st_coordinates(geometry)[, "Y"]
    ) |>
    sf::st_drop_geometry()
}
```

Using insistently we can ensure that temporal failures are tried, this
makes the code slightly more robust for a failed download attempt

``` r
vid_for_timestamp_insistently <- purrr::insistently(
  vid_for_timestamp,
  quiet = FALSE,
  rate = purrr::rate_backoff(
    max_times = 5,
    pause_base = 5,
    jitter = F,
    pause_cap = 60
  )
)
```

``` r
ppi_df <- purrr::map(
  .progress = FALSE,
  animation_timestamps,
  vid_for_timestamp_insistently
)
```

    ## Loading required namespace: raster

``` r
ppi_df_summarized <- ppi_df |>
  bind_rows()
```

Using the `gganimate` package we can create an animation from all
timestamps.

``` r
require(gganimate)
animation <- ggplot() +
  geom_tile(
    data = ppi_df_summarized,
    aes(fill = VID, x = X, y = Y),
    color = NA
  ) +
  scale_fill_viridis_c(
    name = "Bird density [1/km^2]",
    option = "turbo",
    limit = c(2, 800),
    transform = "log10"
  ) +
  geom_sf(data = map, fill = NA, color = "black", size = .3) +
  geom_sf(data = radarLocation) +
  coord_sf(
    expand = F,
    xlim = st_bbox(radarGrid)[c("xmin", "xmax")],
    ylim = st_bbox(radarGrid)[c("ymin", "ymax")]
  ) +
  labs(
    title = "Bird density over the Netherlands",
    subtitle = "Time: {current_frame}",
  ) +
  transition_manual(time)
```

    ## Loading required namespace: raster

This animation can then be controlled and printed in different ways. For
example, as a `gif` or video file.

``` r
animate(
  animation,
  nframes = length(animation_timestamps),
  fps = 5
)
```

    ## Warning in scale_fill_viridis_c(name = "Bird density [1/km^2]", option =
    ## "turbo", : log-10 transformation introduced infinite values.

![](case_study_3_reproducible_files/figure-gfm/printAnimation1R-1.gif)<!-- -->

In this animation we again see the high density of migratory birds along
the coast. This corresponds to observations by a radar (MR1) located
there. At the same time a similar radar located in Amsterdam did not
observe a migration peek, in the animation we see that a high density of
birds migrates south of Amterdam and over the IJselmeer but not throught
Amsterdam.
