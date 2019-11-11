itifeedback <- c(2709,2145,2098,2793,
                 1597,2348,2417,1878,
                 2107,1417,1426,2052,
                 2361,3292,1666,2093,
                 1958,1033,1780,1102,
                 2420,1555,2050,1727,
                 2151,1699,2244,2369,
                 2855,1902,930,1580,
                 2677,1463,2480,2062,
                 2718,1019,1901,1396)

wes_cols= wes_palette("GrandBudapest1", n = 2)
cols <- wes_palette("Darjeeling2", n = 4)

ggplot(data=itifeedback) + aes(itifeedback) + geom_histogram(binwidth=100,col = cols[2], fill=cols[1]) +
scale_x_continuous(name = "Distribution of ITI feedbacks",
                   breaks = seq(750, 3000, 500),
                   limits=c(750, 3000)) +
  scale_y_continuous(name = "Count") +
  geom_rug(alpha=1/2)

itis <- c(2448,2379,2659,2656,
          2067,2484,2417,2813,
          3046,3054,2068,2538,
          1892,1943,2496,3266,
          2115,2685,2387,3058,
          1955,2516,2776,3050,
          3272,2542,1754,2128,
          1969,3675,2192,2874,
          2403,2944,2117,1798,
          1788,2744,2411,2401)