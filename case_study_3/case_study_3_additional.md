Additional exploration
================

<!---# Source - https://stackoverflow.com/a
# Posted by Thomas
# Retrieved 2025-12-09, License - CC BY-SA 3.0 -->

First some code from the reproducible example is repeated to use the
same time frame

``` r
require(bioRad)
require(getRad)
require(lubridate)
require(dplyr)
require(sf)
require(ggplot2)
require(stars)
```

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

Lets take a bit wider interval to see the context

``` r
small_scale_radar_interval<-lubridate::interval(time-hours(10), time+hours(10))
```

## Robin Radars located at sea

To investigate the high densities at sea we analyse two robin radars
located there.

``` r
require(uvaRR)
trks <- purrr::map(
  c("rws01", "borssele"),
  ~ rr_track(
    database = .x,
    ts =c(int_start(small_scale_radar_interval), int_end(small_scale_radar_interval)) ,
    interv = 60 * 60,
    geom = T,
    classification = c("BIRD"),
    movement_metrics = c("dir", "dist", "edist", "spd", "str", "z")
  ) |>
    tibble::add_column(radar = .x)
) |>
  bind_rows() |>
  mutate(trajectory = sf::st_as_sfc(trajectory, crs = 4326)) |>
  st_as_sf() |>
  mutate(
    radarid_unique = trajectory_radarid |>
      purrr::map(
        ~ unlist(purrr::map(
          purrr::map(strsplit(gsub("\\{|\\}", "", .x), ","), as.numeric),
          unique
        ))
      ),
    radarid_unique_chr = purrr::map_chr(
      purrr::map(radarid_unique, as.character),
      paste,
      collapse = " "
    ),
    radar_type = case_match(
      radarid_unique_chr,
      c("11 14", "5 8") ~ "mixed",
      c("14", "8") ~ "vertical",
      c("11", "5") ~ "horizontal"
    )
  )
```

    ## no assignable properties filter specified in WHERE-clause

    ## classification filter included in WHERE-clause

    ## no tracktype filter (RaAz, RaEl, RaAzEl) specified in WHERE-clause

    ## nr_of_plots filter included in WHERE-clause is less than 1500

    ## no polygon filter specified in WHERE-clause

    ## start while loop

    ## 0h 1 min 36 sec

    ## no assignable properties filter specified in WHERE-clause

    ## classification filter included in WHERE-clause

    ## no tracktype filter (RaAz, RaEl, RaAzEl) specified in WHERE-clause

    ## nr_of_plots filter included in WHERE-clause is less than 1500

    ## no polygon filter specified in WHERE-clause

    ## start while loop

    ## 0h 0 min 33 sec

There is an increased density of birds in Luchterduinen (`rws01`) while
this is absent from Borssele.

``` r
ggplot(trks) +
  geom_histogram(aes(x = timestamp_start, fill = classification)) +
  facet_grid(radar ~ .) +
  geom_vline(xintercept = time)
```

    ## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

![](case_study_3_additional_files/figure-gfm/rrplot-1.png)<!-- -->

``` r
ggplot(trks) +
  geom_histogram(aes(x = timestamp_start, fill = radar_type)) +
  facet_grid(radar ~ .) +
  geom_vline(xintercept = time)
```

    ## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

![](case_study_3_additional_files/figure-gfm/rrplot-2.png)<!-- --> A
directional histogram seems to confirm that the movements are vey
directional corresponding to migratory movements

``` r
ggplot(
  trks |>
    dplyr::filter(
      radar_type != "vertical",
      timestamp_start <
        lubridate::ceiling_date(time + lubridate::hours(2), "hours"),
      timestamp_start >
        lubridate::floor_date(time - lubridate::hours(2), "hours")
    )
) +
  geom_histogram(aes(x = dir, fill = radar_type), binwidth = 30) +
  facet_grid(radar ~ lubridate::floor_date(timestamp_start, "hour")) +
  coord_polar() +
  scale_x_continuous(limits = c(0, 360), breaks = 0:3*90)
```

    ## Warning: Removed 32 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](case_study_3_additional_files/figure-gfm/dirhist-1.png)<!-- -->

## MR1 Radars

``` r
bs_radar_databases = c("sbrs_kid_bs291", "sbrs_ams_bs290")
bs_data <- purrr::map(
  bs_radar_databases, ~
    DBI::dbGetQuery(
      as.is = F,
      uvaAuth::get_connection_postgresql(.x),
      glue::glue(
        "select * from collection left join rf_classification on rf_classification.echo = collection.row
      left join rfclasses on rfclasses.id = rf_classification.class  where time_stamp between '{format(lubridate::int_start(small_scale_radar_interval))}' and '{format(lubridate::int_end(small_scale_radar_interval))}';
"
      )
    ) |>
      dplyr::mutate(
        speed = feature3,
        azimuth = feature2,
        altitude_AGL = feature1
      ) %>%
      dplyr::select(-starts_with("feature")) %>%
      tibble::add_column(database = .x)
) |>
  bind_rows()
```

``` r
bs_data_filtered<-bs_data |> dplyr::filter(!name%in%c('insect', 'undefined', 'nonbio','precipitation'))
table(bs_data$name)
```

    ## 
    ##     bird_flock         insect     large_bird         nonbio passerine_type 
    ##           1336           8781             99           1872           9903 
    ##  precipitation     swift_type      undefined      unid_bird     wader_type 
    ##            178            536              1           5295           1623

``` r
ggplot(bs_data_filtered) +
  geom_histogram(aes(x = time_stamp, fill = name)) +
  facet_grid(database ~ .) +
  geom_vline(xintercept = time)
```

    ## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.

![](case_study_3_additional_files/figure-gfm/explore-1.png)<!-- -->

Here we see that the radar in Kijkduin has a migratory peak in the
morning while that is absent from Amsterdam.

Lets explore the directions in Kijkduin:

``` r
ggplot(
  bs_data_filtered |>
    dplyr::filter(
     database=='sbrs_kid_bs291',
      time_stamp <
        lubridate::ceiling_date(time + lubridate::hours(2), "hours"),
      time_stamp >
        lubridate::floor_date(time - lubridate::hours(2), "hours")
    )
) +
  geom_histogram(aes(x = azimuth, fill = name), binwidth = 15) +
  facet_wrap(. ~ lubridate::floor_date(time_stamp, "hour")) +
  coord_polar() +
  scale_x_continuous(limits = c(0, 360), breaks = 0:3*90)
```

    ## Warning: Removed 1487 rows containing non-finite outside the scale range
    ## (`stat_bin()`).

    ## Warning: Removed 44 rows containing missing values or values outside the scale range
    ## (`geom_bar()`).

![](case_study_3_additional_files/figure-gfm/dirhistBs-1.png)<!-- -->
