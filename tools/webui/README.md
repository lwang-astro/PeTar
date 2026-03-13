# PeTar Option Studio

This folder contains a static frontend that turns `petar -h` output into a configurable option UI.

## Files

- `index.html`: main page
- `styles.css`: visual style
- `app.js`: help parser, option renderer, command generator
- `petar-help.txt`: captured output of `petar -h` (used as the default data source)

## Run

From this directory:

```bash
python3 -m http.server 8080
```

Open `http://localhost:8080` in a browser.

## Update options

When PeTar options change, refresh the help source:

```bash
petar -h > petar-help.txt 2>&1
```

Then reload the page.
