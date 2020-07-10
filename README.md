# SunCalc

SunCalc is a tiny ColdFusion (Lucee) library for calculating sun position,
sunlight phases (times for sunrise, sunset, dusk, etc.),
moon position and lunar phase for the given location and time,
created by [Scott Steinbeck](http://scottsteinbeck.com) ([@scottsteinbeck](https://github.com/scottsteinbeck))
as a part of the [SunCalc.net project](http://suncalc.net).

Most calculations are based on the formulas given in the excellent Astronomy Answers articles
about [position of the sun](http://aa.quae.nl/en/reken/zonpositie.html)
and [the planets](http://aa.quae.nl/en/reken/hemelpositie.html).
You can read about different twilight phases calculated by SunCalc
in the [Twilight article on Wikipedia](http://en.wikipedia.org/wiki/Twilight).

## Usage example

```javascript
// get today's sunlight times for London
var times = SunCalc.getTimes(now(), 51.5, -0.1);

// format sunrise time from the Date object
var sunriseStr = times.sunrise.getHours() + ":" + times.sunrise.getMinutes();

// get position of the sun (azimuth and altitude) at today's sunrise
var sunrisePos = SunCalc.getPosition(times.sunrise, 51.5, -0.1);

// get sunrise azimuth in degrees
var sunriseAzimuth = (sunrisePos.azimuth * 180) / Math.PI;
```

SunCalc is also available as an NPM package:

```bash
$ box install suncalc
```

```js
var SunCalc = CreateObject("component", "suncalc");
```

## Reference

### Sunlight times

```javascript
SunCalc.getTimes(
  /*Date*/ date,
  /*Number*/ latitude,
  /*Number*/ longitude,
  /*Number (default=0)*/ height
);
```

Returns an object with the following properties (each is a `Date` object):

| Property        | Description                                                              |
| --------------- | ------------------------------------------------------------------------ |
| `sunrise`       | sunrise (top edge of the sun appears on the horizon)                     |
| `sunriseEnd`    | sunrise ends (bottom edge of the sun touches the horizon)                |
| `goldenHourEnd` | morning golden hour (soft light, best time for photography) ends         |
| `solarNoon`     | solar noon (sun is in the highest position)                              |
| `goldenHour`    | evening golden hour starts                                               |
| `sunsetStart`   | sunset starts (bottom edge of the sun touches the horizon)               |
| `sunset`        | sunset (sun disappears below the horizon, evening civil twilight starts) |
| `dusk`          | dusk (evening nautical twilight starts)                                  |
| `nauticalDusk`  | nautical dusk (evening astronomical twilight starts)                     |
| `night`         | night starts (dark enough for astronomical observations)                 |
| `nadir`         | nadir (darkest moment of the night, sun is in the lowest position)       |
| `nightEnd`      | night ends (morning astronomical twilight starts)                        |
| `nauticalDawn`  | nautical dawn (morning nautical twilight starts)                         |
| `dawn`          | dawn (morning nautical twilight ends, morning civil twilight starts)     |

```javascript
SunCalc.addTime(
  /*Number*/ angleInDegrees,
  /*String*/ morningName,
  /*String*/ eveningName
);
```

Adds a custom time when the sun reaches the given angle to results returned by `SunCalc.getTimes`.

`SunCalc.times` property contains all currently defined times.

### Sun position

```javascript
SunCalc.getPosition(
  /*Date*/ timeAndDate,
  /*Number*/ latitude,
  /*Number*/ longitude
);
```

Returns an object with the following properties:

- `altitude`: sun altitude above the horizon in radians,
  e.g. `0` at the horizon and `PI/2` at the zenith (straight over your head)
- `azimuth`: sun azimuth in radians (direction along the horizon, measured from south to west),
  e.g. `0` is south and `Math.PI * 3/4` is northwest

### Moon position

```javascript
SunCalc.getMoonPosition(
  /*Date*/ timeAndDate,
  /*Number*/ latitude,
  /*Number*/ longitude
);
```

Returns an object with the following properties:

- `altitude`: moon altitude above the horizon in radians
- `azimuth`: moon azimuth in radians
- `distance`: distance to moon in kilometers
- `parallacticAngle`: parallactic angle of the moon in radians

### Moon illumination

```javascript
SunCalc.getMoonIllumination(/*Date*/ timeAndDate);
```

Returns an object with the following properties:

- `fraction`: illuminated fraction of the moon; varies from `0.0` (new moon) to `1.0` (full moon)
- `phase`: moon phase; varies from `0.0` to `1.0`, described below
- `angle`: midpoint angle in radians of the illuminated limb of the moon reckoned eastward from the north point of the disk;
  the moon is waxing if the angle is negative, and waning if positive

Moon phase value should be interpreted like this:

| Phase | Name            |
| ----: | --------------- |
|     0 | New Moon        |
|       | Waxing Crescent |
|  0.25 | First Quarter   |
|       | Waxing Gibbous  |
|   0.5 | Full Moon       |
|       | Waning Gibbous  |
|  0.75 | Last Quarter    |
|       | Waning Crescent |

By subtracting the `parallacticAngle` from the `angle` one can get the zenith angle of the moons bright limb (anticlockwise).
The zenith angle can be used do draw the moon shape from the observers perspective (e.g. moon lying on its back).

### Moon rise and set times

```js
SunCalc.getMoonTimes(/*Date*/ date, /*Number*/ latitude, /*Number*/ longitude[, inUTC])
```

Returns an object with the following properties:

- `rise`: moonrise time as `Date`
- `set`: moonset time as `Date`
- `alwaysUp`: `true` if the moon never rises/sets and is always _above_ the horizon during the day
- `alwaysDown`: `true` if the moon is always _below_ the horizon

By default, it will search for moon rise and set during local user's day (frou 0 to 24 hours).
If `inUTC` is set to true, it will instead search the specified date from 0 to 24 UTC hours.

## Changelog

#### 1.0.0 &mdash; Jul 10, 2020

- First commit.
