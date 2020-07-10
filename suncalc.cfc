component {

	variables.Math = createObject("java","java.lang.Math");
	variables.PI   = Math.PI;
	variables.rad  = PI / 180;
	private numeric function atan ( y, x ){
		return Math.atan2( javacast('double', y), javacast('double', x))
	}

	variables.J0 = 0.0009;
	// date/time constants and conversions
	variables.dayMs = 1000 * 60 * 60 * 24;
	variables.J1970 = 2440588;
	variables.J2000 = 2451545;
	variables.e = rad * 23.4397; // obliquity of the Earth

	// sun times configuration (angle, morning name, evening name)

	variables.times = [
		[-0.833, 'sunrise',       'sunset'      ],
		[  -0.3, 'sunriseEnd',    'sunsetStart' ],
		[    -6, 'dawn',          'dusk'        ],
		[   -12, 'nauticalDawn',  'nauticalDusk'],
		[   -18, 'nightEnd',      'night'       ],
		[     6, 'goldenHourEnd', 'goldenHour'  ]
	];

	// adds a custom time to the times config

	public void function addTime (angle, riseName, setName) {
		variables.times.append([angle, riseName, setName]);
	};

	// calculates sun position for a given date and latitude/longitude

	public struct function getPosition (date, lat, lng) {

		var lw  = rad * -lng;
		var phi = rad * lat;
		var d   = toDays(date);

		var c  = sunCoords(d);
		var H  = siderealTime(d, lw) - c.ra;

		return {
			'azimuth': azimuth(H, phi, c.dec),
			'altitude': altitude(H, phi, c.dec)
		};
	};



	// calculates sun times for a given date, latitude/longitude, and, optionally,
	// the observer height (in meters) relative to the horizon
	public struct function getTimes (date, lat, lng, height = 0) {

		var lw = rad * -lng;
		var phi = rad * lat;

		var dh = observerAngle(height);

		var d = toDays(date);
		var n = julianCycle(d, lw);
		var ds = approxTransit(0, lw, n);

		var M = solarMeanAnomaly(ds);
		var L = eclipticLongitude(M);
		var dec = declination(L, 0);

		var Jnoon = solarTransitJ(ds, M, L);

		var result = {
			'solarNoon': fromJulian(Jnoon),
			'nadir': fromJulian(Jnoon - 0.5)
		};

		var  len = times.len();
		for (var i = 1; i <= len; i += 1) {
			var time = times[i];
			var h0 = (time[1] + dh) * rad;

			var Jset = getSetJ(h0, lw, phi, dec, n, M, L);
			var Jrise = Jnoon - (Jset - Jnoon);

			result[time[2]] = fromJulian(Jrise);
			result[time[3]] = fromJulian(Jset);
		}

		return result;
	};



	public struct function getMoonPosition (date, lat, lng) {

		var lw  = rad * -lng;
		var phi = rad * lat;
		var d   = toDays(date);

		var c = moonCoords(d);
		var H = siderealTime(d, lw) - c.ra;
		var h = altitude(H, phi, c.dec);
			// formula 14.1 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
		var pa = atan(sin(H), tan(phi) * cos(c.dec) - sin(c.dec) * cos(H));

		var h = h + astroRefraction(h); // altitude correction for refraction

		return {
			'azimuth': azimuth(H, phi, c.dec),
			'altitude': h,
			'distance': c.dist,
			'parallacticAngle': pa
		};
	};


	// calculations for illumination parameters of the moon,
	// based on http://idlastro.gsfc.nasa.gov/ftp/pro/astro/mphase.pro formulas and
	// Chapter 48 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
	public struct function getMoonIllumination (date = now()) {

		var d = toDays(date);
		var s = sunCoords(d);
		var m = moonCoords(d);

		var sdist = 149598000; // distance from Earth to Sun in km

		var phi = acos(sin(s.dec) * sin(m.dec) + cos(s.dec) * cos(m.dec) * cos(s.ra - m.ra));
		var inc = atan(sdist * sin(phi), m.dist - sdist * cos(phi));
		var angle = atan(cos(s.dec) * sin(s.ra - m.ra), sin(s.dec) * cos(m.dec) -
					cos(s.dec) * sin(m.dec) * cos(s.ra - m.ra));

		return {
			'fraction': (1 + cos(inc)) / 2,
			'phase': 0.5 + 0.5 * inc * (angle < 0 ? -1 : 1) / Math.PI,
			'angle': angle
		};
	};


	private date function hoursLater(date, h) {
		return createObject("java", "java.util.Date").init(javaCast("long", date.getTime() + h * dayMs / 24));
	}


	// calculations for moon rise/set times are based on http://www.stargazing.net/kepler/moonrise.html article
	public struct function getMoonTimes (date = now(), lat, lng) {
		var t = createDate(year(date),month(date),day(date));
		var hc = 0.133 * rad;
		var h0 = getMoonPosition(t, lat, lng).altitude - hc;

		// go in 2-hour chunks, each time seeing if a 3-point quadratic curve crosses zero (which means rise or set)
		for (var i = 1; i <= 24; i += 2) {
		var h1 = getMoonPosition(hoursLater(t, i), lat, lng).altitude - hc;
		var h2 = getMoonPosition(hoursLater(t, i + 1), lat, lng).altitude - hc;

		var a = (h0 + h2) / 2 - h1;
		var b = (h2 - h0) / 2;
		var xe = -b / (2 * a);
		var ye = (a * xe + b) * xe + h1;
		var d = b * b - 4 * a * h1;
		var roots = 0;
		var rise = 0;
		var set = 0;

			if (d >= 0) {
				dx = Math.sqrt(d) / (Math.abs(a) * 2);
				x1 = xe - dx;
				x2 = xe + dx;
				if (Math.abs(x1) <= 1) roots++;
				if (Math.abs(x2) <= 1) roots++;
				if (x1 < -1) x1 = x2;
			}

			if (roots === 1) {
				if (h0 < 0) rise = i + x1;
				else set = i + x1;

			} else if (roots === 2) {
				rise = i + (ye < 0 ? x2 : x1);
				set = i + (ye < 0 ? x1 : x2);
			}

			if (rise && set) break;

			h0 = h2;
		}

		var result = {};

		if (rise) result.rise = hoursLater(t, rise);
		if (set) result.set = hoursLater(t, set);

		if (!rise && !set) result[ye > 0 ? 'alwaysUp' : 'alwaysDown'] = true;

		return result;
	};


	// moon calculations, based on http://aa.quae.nl/en/reken/hemelpositie.html formulas
	public struct function moonCoords(d) { // geocentric ecliptic coordinates of the moon

		var L = rad * (218.316 + 13.176396 * d); // ecliptic longitude
		var M = rad * (134.963 + 13.064993 * d); // mean anomaly
		var F = rad * (93.272 + 13.229350 * d);  // mean distance

		var l  = L + rad * 6.289 * sin(M); // longitude
		var b  = rad * 5.128 * sin(F);     // latitude
		var dt = 385001 - 20905 * cos(M);  // distance to the moon in km

		return {
			'ra': rightAscension(l, b),
			'dec': declination(l, b),
			'dist': dt
		};
	}



	// sun calculations are based on http://aa.quae.nl/en/reken/zonpositie.html formulas
	private numeric function toJulian(date) { return date.getTime() / dayMs - 0.5 + J1970; }
	private numeric function fromJulian(j)  { return createObject("java", "java.util.Date").init(javaCast("long", ((j + 0.5 - J1970) * dayMs))); }
	private numeric function toDays(date)   { return toJulian(date) - J2000; }


	// general calculations for position
	private numeric function rightAscension(l, b) { return atan(sin(l) * cos(e) - tan(b) * sin(e), cos(l)); }
	private numeric function declination(l, b)    { return asin(sin(b) * cos(e) + cos(b) * sin(e) * sin(l)); }

	private numeric function azimuth(H, phi, dec)  { return atan(sin(H), cos(H) * sin(phi) - tan(dec) * cos(phi)); }
	private numeric function altitude(H, phi, dec) { return asin(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(H)); }

	private numeric function siderealTime(d, lw) { return rad * (280.16 + 360.9856235 * d) - lw; }

	private numeric function astroRefraction(h) {
		if (h < 0) // the following formula works for positive altitudes only.
			h = 0; // if h = -0.08901179 a div/0 would occur.

		// formula 16.4 of "Astronomical Algorithms" 2nd edition by Jean Meeus (Willmann-Bell, Richmond) 1998.
		// 1.02 / tan(h + 10.26 / (h + 5.10)) h in degrees, result in arc minutes -> converted to rad:
		return 0.0002967 / Math.tan(h + 0.00312536 / (h + 0.08901179));
	}

	// general sun calculations
	private numeric function solarMeanAnomaly(d) { return rad * (357.5291 + 0.98560028 * d); }
	private numeric function eclipticLongitude(M) {

		var C = rad * (1.9148 * sin(M) + 0.02 * sin(2 * M) + 0.0003 * sin(3 * M)); // equation of center
			P = rad * 102.9372; // perihelion of the Earth

		return M + C + P + PI;
	}

	private struct function sunCoords(d) {

		var M = solarMeanAnomaly(d);
		var L = eclipticLongitude(M);

		return {
			dec: declination(L, 0),
			ra: rightAscension(L, 0)
		};
	}

	// calculations for sun times

	private numeric function julianCycle(d, lw) { return Math.round(d - J0 - lw / (2 * PI)); }

	private numeric function approxTransit(Ht, lw, n) { return J0 + (Ht + lw) / (2 * PI) + n; }
	private numeric function solarTransitJ(ds, M, L)  { return J2000 + ds + 0.0053 * sin(M) - 0.0069 * sin(2 * L); }

	private numeric function hourAngle(h, phi, d) { return acos((sin(h) - sin(phi) * sin(d)) / (cos(phi) * cos(d))); }
	private numeric function observerAngle(height) { return -2.076 * Math.sqrt(height) / 60; }

	// returns set time for the given sun altitude
	private numeric function getSetJ(h, lw, phi, dec, n, M, L) {

		var w = hourAngle(h, phi, dec);
		var a = approxTransit(w, lw, n);
		return solarTransitJ(a, M, L);
	}
}