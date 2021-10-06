qopen go --filter-events '["magnitude > 2.5"]' --prefix "01_go/"
qopen fixed --input "01_go/results.json" --align-sites --prefix "02_sites/"
qopen source --input "01_go/results.json" --input-sites "02_sites/results.json" --prefix "03_mag/" --overwrite-conf '{"skip": {"coda_window": 1, "num_pairs": 2}}'
qopen recalc_source --input "03_mag/results.json" --prefix "04_mag_nconst/" --seismic-moment-options '{"fc": null, "n": 1.88, "gamma": 2, "num_points": 4}'
qopen go -e 20186784 --dump-fitpkl 01_go/fits_%s.pkl --no-plots --output null
qopen go -e 20186784 --prefix 01_go_20186784/ --plot-optimization --overwrite-conf '{"g0_bounds": [1e-6, 1e-4]}'
