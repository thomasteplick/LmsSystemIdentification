/*
    Create an LMS System Identifier or System Modeler.
	Plot the impulse and frequency responses of the adaptive
	filter or Plant.

	The html/template package is used to generate the html sent to the client.
	Use CSS display grid to display a 300x300 grid of cells.
	Use CSS flexbox to display the labels on the x and y axes.

	Calculate the power spectral density (PSD) using Welch periodogram spectral
	estimation and plot it.  Average the periodograms using 50% overlap with a
	user-chosen window function such as Bartlett, Hanning, Hamming, or Welch.
*/

package main

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"math/cmplx"
	"math/rand"
	"net/http"
	"os"
	"path"
	"strconv"
	"strings"
	"sync"
	"text/template"

	"github.com/mjibson/go-dsp/fft"
)

const (
	rows                       = 300                                  // #rows in grid
	columns                    = 300                                  // #columns in grid
	block                      = 512                                  // size of buf1 and buf2, chunks of data to process
	tmpltime                   = "templates/plottimedata.html"        // html template address
	tmplfrequency              = "templates/plotfrequencydata.html"   // html template address
	tmplplotresponse           = "templates/plotresponse.html"        // html template address
	tmplLmsSystemIdentifier    = "templates/LmsSystemIdentifier.html" // html template address
	tmplfiltersignal           = "templates/filtersignal.html"        // html template address
	addr                       = "127.0.0.1:8080"                     // http server listen address
	patternFilterSignal        = "/filtersignal"                      // http handler pattern for filtering using the LMS System Identifier
	patternLmsSystemIdentifier = "/lmssystemidentifier"               // http handler pattern for creating LMS System Identifier
	patternPlotResponse        = "/plotresponse"                      // http handler for plotting impulse or frequency responses
	xlabels                    = 11                                   // # labels on x axis
	ylabels                    = 11                                   // # labels on y axis
	dataDir                    = "data/"                              // directory for the data files
	deg2rad                    = math.Pi / 180.0                      // convert degrees to radians
	signal                     = "signal.txt"                         // signal to plot
	lmsSystemIdentifierFile    = "lmssystemidentifier.txt"            // LMS adaptive filter for System Identification
	plantModelFile             = "plantmodel.txt"                     // plant difference equation coefficients
	plantImpulseResponseFile   = "plantImpulseResponse.txt"           // plant Impulse Response
	twoPi                      = 2.0 * math.Pi                        // Two PI ~ 6.28
)

// Type to contain all the HTML template actions
type PlotT struct {
	Grid        []string // plotting grid
	Status      string   // status of the plot
	Xlabel      []string // x-axis labels
	Ylabel      []string // y-axis labels
	Filename    string   // filename to plot
	SampleFreq  string   // data sampling rate in Hz
	FFTSegments string   // FFT segments, K
	FFTSize     string   // FFT size
	Samples     string   // complex samples in data file
}

// Type to hold the minimum and maximum data values
type Endpoints struct {
	xmin float64
	xmax float64
	ymin float64
	ymax float64
}

// Window function type
type Window func(n int, m int) complex128

// data to be supplied between the pipeline stages through the channel
type ComChan struct {
	desired float64 // primary source signal to the LMS algorithm
	in      float64 // input to the stage in the pipeline
}

// LMS algorithm  data to create the system identifier
type LMSAlgorithm struct {
	gain       float64   // gain of filter
	trials     int       // number of trials for algorithm
	order      int       // adaptive filter order
	samples    int       // number of samples
	samplerate int       // sample rate in Hz
	wEnsemble  []float64 // ensemble average coefficients
	wTrial     []float64 // trial coefficients
	wg         sync.WaitGroup
	toChan     chan float64 // synchronized channel to plant from signal source generator
	fromChan   chan ComChan // synchronized channel to adaptive filter from plant
	poleRad    []float64    // plant pole radius ->  r*cos(ang) + j r*sin(ang)
	poleAng    []float64    // plant pole angle in degrees
	zeroAng    []float64    // plant zero angle in degrees
	snr        int          // signal-to-noise ratio for noisy channel
}

// biquads are second-order zero-pole sections consisting of a
// numerator and denominator polynomial in z^(-1)
type Biquad struct {
	num   [3]float64
	denom [3]float64
}

var (
	timeTmpl                *template.Template
	freqTmpl                *template.Template
	lmsSystemIdentifierTmpl *template.Template
	plotresponseTmpl        *template.Template
	winType                 = []string{"Bartlett", "Welch", "Hamming", "Hanning"}
)

// Bartlett window
func bartlett(n int, m int) complex128 {
	real := 1.0 - math.Abs((float64(n)-float64(m))/float64(m))
	return complex(real, 0)
}

// Welch window
func welch(n int, m int) complex128 {
	x := math.Abs((float64(n) - float64(m)) / float64(m))
	real := 1.0 - x*x
	return complex(real, 0)
}

// Hamming window
func hamming(n int, m int) complex128 {
	return complex(.54-.46*math.Cos(math.Pi*float64(n)/float64(m)), 0)
}

// Hanning window
func hanning(n int, m int) complex128 {
	return complex(.5-.5*math.Cos(math.Pi*float64(n)/float64(m)), 0)
}

// Rectangle window
func rectangle(n int, m int) complex128 {
	return 1.0
}

// init parses the html template files and is done only once at startup
func init() {
	timeTmpl = template.Must(template.ParseFiles(tmpltime))
	freqTmpl = template.Must(template.ParseFiles(tmplfrequency))
	plotresponseTmpl = template.Must(template.ParseFiles(tmplplotresponse))
	lmsSystemIdentifierTmpl = template.Must(template.ParseFiles(tmplLmsSystemIdentifier))
}

// findEndpointsMA finds the minimum and maximum data values
func (ep *Endpoints) findEndpointsMA(input *bufio.Scanner) {
	ep.xmax = -math.MaxFloat64
	ep.xmin = math.MaxFloat64
	ep.ymax = -math.MaxFloat64
	ep.ymin = math.MaxFloat64

	n := 0

	for input.Scan() {
		// Each line has 1 value
		line := input.Text()
		var (
			y   float64
			err error
		)
		// MA impulse response plot
		if y, err = strconv.ParseFloat(line, 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", line, err)
			continue
		}
		n++

		if y > ep.ymax {
			ep.ymax = y
		}
		if y < ep.ymin {
			ep.ymin = y
		}
	}
	// impulse response plot using sample number for the x-axis
	ep.xmin = 0.0
	ep.xmax = float64(n - 1)
}

// findEndpointsARMA finds the minimum and maximum data values
func (ep *Endpoints) findEndpointsARMA(input *bufio.Scanner) {
	ep.xmax = -math.MaxFloat64
	ep.xmin = math.MaxFloat64
	ep.ymax = -math.MaxFloat64
	ep.ymin = math.MaxFloat64

	var (
		err  error
		temp float64
	)
	a := make([]float64, 0)
	b := make([]float64, 0)

	// Save the impulse response
	impulsefile, err := os.Create(path.Join(dataDir, plantImpulseResponseFile))
	if err != nil {
		fmt.Printf("Create %s error: %v", plantImpulseResponseFile, err)
		return
	}
	defer impulsefile.Close()

	// Retrieve the plant model
	for input.Scan() {
		// Each line has 2 values
		line := input.Text()
		values := strings.Split(line, " ")
		// ARMA impulse response plot
		if temp, err = strconv.ParseFloat(values[0], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
			continue
		}
		b = append(b, temp)
		if temp, err = strconv.ParseFloat(values[1], 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
			continue
		}
		a = append(a, temp)
	}
	ncoeff := len(a)

	x := make([]float64, ncoeff)
	y := make([]float64, ncoeff)

	for i := 0; i < ncoeff; i++ {
		temp = 0.0
		x[i] = 1
		k := (i - 1) % ncoeff
		if k < 0 {
			k = ncoeff + k
		}
		x[k] = 0
		for j := 0; j < ncoeff; j++ {
			k = (i - j) % ncoeff
			if k < 0 {
				k = ncoeff + k
			}
			temp += b[j] * x[k]
		}
		for j := 1; j < ncoeff; j++ {
			k = (i - j) % ncoeff
			if k < 0 {
				k = ncoeff + k
			}
			temp -= a[j] * y[k]
		}
		y[i] = temp

		// Save to plantImpulseResponseFile
		fmt.Fprintf(impulsefile, "%f\n", y[i])

		if temp > ep.ymax {
			ep.ymax = temp
		}
		if temp < ep.ymin {
			ep.ymin = temp
		}
	}

	i := 0
	const its int = 100
	for n := 0; n < its; n++ {
		temp = 0.0
		for j := 1; j < ncoeff; j++ {
			k := (i - j) % ncoeff
			if k < 0 {
				k = ncoeff + k
			}
			temp -= a[j] * y[k]
		}
		y[i] = temp

		// Save to plantImpulseResponseFile
		fmt.Fprintf(impulsefile, "%f\n", y[i])

		i = (i + 1) % ncoeff

		if temp > ep.ymax {
			ep.ymax = temp
		}
		if temp < ep.ymin {
			ep.ymin = temp
		}
	}

	// impulse response plot using sample number as the x-axis
	ep.xmin = 0.0
	ep.xmax = float64(ncoeff + its - 1)
}

// gridFill inserts the data points in the grid
func gridFill(plot *PlotT, xscale float64, yscale float64, endpoints Endpoints, input *bufio.Scanner) error {
	var x float64 = -1
	for input.Scan() {
		line := input.Text()
		// Each line has 1 value
		var (
			y   float64
			err error
		)
		x++
		if y, err = strconv.ParseFloat(line, 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", line, err)
			return err
		}

		// Check if inside the zoom values
		if x < endpoints.xmin || x > endpoints.xmax || y < endpoints.ymin || y > endpoints.ymax {
			continue
		}

		// This cell location (row,col) is on the line
		row := int((endpoints.ymax-y)*yscale + .5)
		col := int((x-endpoints.xmin)*xscale + .5)
		plot.Grid[row*columns+col] = "online"
	}
	return nil
}

// gridFillInterp inserts the data points in the grid and draws a straight line between points
func gridFillInterp(plot *PlotT, xscale float64, yscale float64, endpoints Endpoints, input *bufio.Scanner) error {

	var (
		x, y         float64
		prevX, prevY float64
		prevSet      bool = true
		err          error
	)

	const lessen = 1
	const increase = 10

	// Get first sample
	input.Scan()
	line := input.Text()
	// Each line has 1 values
	x = 0
	if y, err = strconv.ParseFloat(line, 64); err != nil {
		fmt.Printf("String %s conversion to float error: %v\n", line, err)
		return err
	}

	// Check if inside the zoom values
	if x < endpoints.xmin || x > endpoints.xmax || y < endpoints.ymin || y > endpoints.ymax {
		prevSet = false
	} else {
		// This cell location (row,col) is on the line
		row := int((endpoints.ymax-y)*yscale + .5)
		col := int((x-endpoints.xmin)*xscale + .5)
		plot.Grid[row*columns+col] = "online"

		prevX = x
		prevY = y
	}

	// Scale factor to determine the number of interpolation points
	lenEP := math.Sqrt((endpoints.xmax-endpoints.xmin)*(endpoints.xmax-endpoints.xmin) +
		(endpoints.ymax-endpoints.ymin)*(endpoints.ymax-endpoints.ymin))

	// Continue with the rest of the points in the file
	for input.Scan() {
		line = input.Text()
		// Each line has 1 value
		x++
		if y, err = strconv.ParseFloat(line, 64); err != nil {
			fmt.Printf("String %s conversion to float error: %v\n", line, err)
			return err
		}

		// Check if inside the zoom values
		if x < endpoints.xmin || x > endpoints.xmax || y < endpoints.ymin || y > endpoints.ymax {
			continue
		} else if !prevSet {
			prevSet = true
			prevX = x
			prevY = y
		}

		// This cell location (row,col) is on the line
		row := int((endpoints.ymax-y)*yscale + .5)
		col := int((x-endpoints.xmin)*xscale + .5)
		plot.Grid[row*columns+col] = "online"

		// Interpolate the points between previous point and current point

		lenEdge := math.Sqrt((x-prevX)*(x-prevX) + (y-prevY)*(y-prevY))
		ncells := increase * int(columns*lenEdge/lenEP) / lessen // number of points to interpolate
		stepX := (x - prevX) / float64(ncells)
		stepY := (y - prevY) / float64(ncells)

		// loop to draw the points
		interpX := prevX
		interpY := prevY
		for i := 0; i < ncells; i++ {
			row := int((endpoints.ymax-interpY)*yscale + .5)
			col := int((interpX-endpoints.xmin)*xscale + .5)
			plot.Grid[row*columns+col] = "online"
			interpX += stepX
			interpY += stepY
		}

		// Update the previous point with the current point
		prevX = x
		prevY = y
	}
	return nil
}

// processTimeDomain plots the time domain data from disk file
func processTimeDomain(w http.ResponseWriter, r *http.Request, filename string) error {

	// main data structure
	var (
		plot      PlotT
		xscale    float64
		yscale    float64
		endpoints Endpoints
	)

	plot.Grid = make([]string, rows*columns)
	plot.Xlabel = make([]string, xlabels)
	plot.Ylabel = make([]string, ylabels)

	// Open file
	f, err := os.Open(filename)
	if err == nil {
		// Mark the data x-y coordinate online at the corresponding
		// grid row/column.
		input := bufio.NewScanner(f)

		// Determine if MA or ARMA, the adaptive filter is MA
		// the plant model is ARMA
		if filename == lmsSystemIdentifierFile {
			endpoints.findEndpointsMA(input)
		} else if filename == plantModelFile {
			endpoints.findEndpointsARMA(input)
		} else {
			f.Close()
			return fmt.Errorf("%s is not a MA/ARMA filter coefficient file", filename)
		}

		f.Close()
		if filename == plantModelFile {
			filename = plantImpulseResponseFile
		}
		f, err = os.Open(filename)
		if err == nil {
			defer f.Close()
			input := bufio.NewScanner(f)

			// Determine if zoom requested and validate endpoints
			zoomxstart := r.FormValue("zoomxstart")
			zoomxend := r.FormValue("zoomxend")
			zoomystart := r.FormValue("zoomystart")
			zoomyend := r.FormValue("zoomyend")
			if len(zoomxstart) > 0 && len(zoomxend) > 0 &&
				len(zoomystart) > 0 && len(zoomyend) > 0 {
				x1, err1 := strconv.ParseFloat(zoomxstart, 64)
				x2, err2 := strconv.ParseFloat(zoomxend, 64)
				y1, err3 := strconv.ParseFloat(zoomystart, 64)
				y2, err4 := strconv.ParseFloat(zoomyend, 64)

				if err1 != nil || err2 != nil || err3 != nil || err4 != nil {
					plot.Status = "Zoom x or y values are not numbers."
					fmt.Printf("Zoom error: x start error = %v, x end error = %v\n", err1, err2)
					fmt.Printf("Zoom error: y start error = %v, y end error = %v\n", err3, err4)
				} else {
					if (x1 < endpoints.xmin || x1 > endpoints.xmax) ||
						(x2 < endpoints.xmin || x2 > endpoints.xmax) || (x1 >= x2) {
						plot.Status = "Zoom values are not in x range."
						fmt.Printf("Zoom error: start or end value not in x range.\n")
					} else if (y1 < endpoints.ymin || y1 > endpoints.ymax) ||
						(y2 < endpoints.ymin || y2 > endpoints.ymax) || (y1 >= y2) {
						plot.Status = "Zoom values are not in y range."
						fmt.Printf("Zoom error: start or end value not in y range.\n")
						fmt.Printf("y1=%v, y2=%v, endpoints.ymin=%v, endpoints.ymax=%v\n",
							y1, y2, endpoints.ymin, endpoints.ymax)
					} else {
						// Valid Zoom endpoints, replace the previous min and max values
						endpoints.xmin = x1
						endpoints.xmax = x2
						endpoints.ymin = y1
						endpoints.ymax = y2
					}
				}
			}

			// Calculate scale factors for x and y
			xscale = (columns - 1) / (endpoints.xmax - endpoints.xmin)
			yscale = (rows - 1) / (endpoints.ymax - endpoints.ymin)

			// Check for interpolation and fill in the grid with the data points
			interp := r.FormValue("interpolate")
			if interp == "interpolate" {
				err = gridFillInterp(&plot, xscale, yscale, endpoints, input)
			} else {
				err = gridFill(&plot, xscale, yscale, endpoints, input)
			}
			if err != nil {
				return err
			}

			// Set plot status if no errors
			if len(plot.Status) == 0 {
				plot.Status = fmt.Sprintf("Status: Data plotted from (%.3f,%.3f) to (%.3f,%.3f)",
					endpoints.xmin, endpoints.ymin, endpoints.xmax, endpoints.ymax)
			}

		} else {
			// Set plot status
			fmt.Printf("Error opening file %s: %v\n", filename, err)
			return fmt.Errorf("error opening file %s: %v", filename, err)
		}
	} else {
		// Set plot status
		fmt.Printf("Error opening file %s: %v\n", filename, err)
		return fmt.Errorf("error opening file %s: %v", filename, err)
	}

	// Construct x-axis labels
	incr := (endpoints.xmax - endpoints.xmin) / (xlabels - 1)
	x := endpoints.xmin
	// First label is empty for alignment purposes
	for i := range plot.Xlabel {
		plot.Xlabel[i] = fmt.Sprintf("%.2f", x)
		x += incr
	}

	// Construct the y-axis labels
	incr = (endpoints.ymax - endpoints.ymin) / (ylabels - 1)
	y := endpoints.ymin
	for i := range plot.Ylabel {
		plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
		y += incr
	}

	// Enter the filename in the form
	plot.Filename = path.Base(filename)

	// Write to HTTP using template and grid
	if err := timeTmpl.Execute(w, plot); err != nil {
		log.Fatalf("Write to HTTP output using template with error: %v\n", err)
	}

	return nil
}

// processFrequencyDomain calculates the power spectral density (PSD) and plots it
func processFrequencyDomain(w http.ResponseWriter, r *http.Request, filename string) error {
	// Use complex128 for FFT computation
	// Get the number of complex samples nn, open file and count lines, close the file

	var (
		plot          PlotT // main data structure to execute with parsed html template
		endpoints     Endpoints
		N             int                                                        //  complex FFT size
		nn            int                                                        // number of complex samples in the data file
		K             int                                                        //  number of segments used in PSD with 50% overlap
		m             int                                                        // complex segment size
		win           string                                                     // FFT window type
		window        = make(map[string]Window, len(winType))                    // map of window functions
		sumWindow     float64                                                    // sum of squared window values for normalization
		normalizerPSD float64                                                    // normalizer for PSD
		PSD           []float64                                                  // power spectral density
		psdMax        float64                                 = -math.MaxFloat64 // maximum PSD value
		psdMin        float64                                 = math.MaxFloat64  // minimum PSD value
		xscale        float64                                                    // data to grid in x direction
		yscale        float64                                                    // data to grid in y direction
		samplingRate  float64                                                    // sampling rate in Hz
		ARMA          bool                                                       // ARMA model for the system
	)

	plot.Grid = make([]string, rows*columns)
	plot.Xlabel = make([]string, xlabels)
	plot.Ylabel = make([]string, ylabels)

	// Put the window functions in the map
	window["Bartlett"] = bartlett
	window["Welch"] = welch
	window["Hamming"] = hamming
	window["Hanning"] = hanning
	window["Rectangle"] = rectangle

	// Need plant model for frequency response, not impulse response
	if filename == plantImpulseResponseFile {
		filename = plantModelFile
	}

	// Open file
	f, err := os.Open(filename)
	if err == nil {
		input := bufio.NewScanner(f)
		// Number of real coefficients
		for input.Scan() {
			line := input.Text()
			if len(line) > 0 {
				nn++
			}
		}
		f.Close()
		fmt.Printf("Data file %s has %d samples\n", filename, nn)
		// make even number of samples so if segments = 1, we won't
		// do the last FFT with one sample
		if nn%2 == 1 {
			nn++
		}

		// Get number of segments from HTML form
		// Number of segments to average the periodograms to reduce the variance
		tmp := r.FormValue("fftsegments")
		if len(tmp) == 0 {
			return fmt.Errorf("enter number of FFT segments")
		}
		K, err = strconv.Atoi(tmp)
		if err != nil {
			fmt.Printf("FFT segments string convert error: %v\n", err)
			return fmt.Errorf("fft segments string convert error: %v", err)
		}

		// Require 1 <= K <= 20
		if K < 1 {
			K = 1
		} else if K > 20 {
			K = 20
		}

		// segment size complex samples
		m = nn / (K + 1)

		// Get window type:  Bartlett, Welch, Hanning, Hamming, etc
		// Multiply the samples by the window to reduce spectral leakage
		// caused by high sidelobes in rectangular window
		win = r.FormValue("fftwindow")
		if len(win) == 0 {
			return fmt.Errorf("enter FFT window type")
		}
		w, ok := window[win]
		if !ok {
			fmt.Printf("Invalid FFT window type: %v\n", win)
			return fmt.Errorf("invalid FFT window type: %v", win)
		}
		// sum the window values for PSD normalization due to windowing
		for i := 0; i < 2*m; i++ {
			x := cmplx.Abs(w(i, m))
			sumWindow += x * x
		}
		fmt.Printf("%s window sum = %.2f\n", win, sumWindow)

		// Get FFT size from HTML form
		// Check FFT Size >= 2*m, using 50%  overlap of segments
		// Check FFT Size is a power of 2:  2^n
		tmp = r.FormValue("fftsize")
		if len(tmp) == 0 {
			return fmt.Errorf("enter FFT size")
		}
		N, err = strconv.Atoi(tmp)
		if err != nil {
			return fmt.Errorf("fft size string convert error: %v", err)
		}

		if N < rows {
			fmt.Printf("FFT size < %d\n", rows)
			N = rows
		} else if N > rows*rows {
			fmt.Printf("FFT size > %d\n", rows*rows)
			N = rows * rows
		}
		// This rounds up to nearest FFT size that is a power of 2
		N = int(math.Exp2(float64(int(math.Log2(float64(N)) + .5))))
		fmt.Printf("N=%v\n", N)

		if N < 2*m {
			fmt.Printf("FFT Size %d not greater than 2%d\n", N, 2*m)
			return fmt.Errorf("fft Size %d not greater than 2*%d", N, 2*m)
		}

		// Power Spectral Density, PSD[N/2] is the Nyquist critical frequency
		// It is the (sampling frequency)/2, the highest non-aliased frequency
		PSD = make([]float64, N/2+1)

		// Reopen the coefficients file
		f, err = os.Open(filename)
		if err == nil {
			defer f.Close()
			bufm_a := make([]complex128, m) // denominator in ARMA
			bufm_b := make([]complex128, m) // numerator in ARMA or MA
			bufN_a := make([]complex128, N) // denominator in ARMA
			bufN_b := make([]complex128, N) // numerator in ARMA or MA
			input := bufio.NewScanner(f)
			// Read in initial m samples into buf[m] to start the processing loop
			diff := 0.0
			for k := 0; k < m; k++ {
				input.Scan()
				line := input.Text()
				// Each line has 1 or 2 space-separated values
				values := strings.Split(line, " ")
				// adaptive filter coefficients MA model, FIR filter
				if len(values) == 1 {
					ARMA = false
					var b float64
					if b, err = strconv.ParseFloat(values[0], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
						continue
					}
					if k != 0 {
						diff += .5
					}
					bufm_b[k] = complex(b, 0)
				} else if len(values) == 2 {
					ARMA = true
					// plant model a and b coefficients, ARMA model
					var a, b float64

					if b, err = strconv.ParseFloat(values[0], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
						continue
					}
					if a, err = strconv.ParseFloat(values[1], 64); err != nil {
						fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
						continue
					}

					if k != 0 {
						diff += .5
					}
					bufm_a[k] = complex(a, 0)
					bufm_b[k] = complex(b, 0)
				}
			}
			// Average the time steps and invert to get the sampling rate
			samplingRate = 1.0 / (diff / float64(m-1))
			fmt.Printf("sampling rate = %.2f\n", samplingRate)

			scanOK := true
			// loop over the rest of the file, reading in m samples at a time until EOF
			for {
				// Put the previous m samples in the front of the buffer N to
				// overlap segments
				if ARMA {
					copy(bufN_a, bufm_a)
				}
				copy(bufN_b, bufm_b)
				// Put the next m samples in back of the previous m samples
				kk := 0
				for k := 0; k < m; k++ {
					scanOK = input.Scan()
					if !scanOK {
						break
					}
					line := input.Text()
					// Each line has 1 - 2 values, depending on MA or ARMA model
					values := strings.Split(line, " ")
					if len(values) == 1 {
						var b float64
						if b, err = strconv.ParseFloat(values[0], 64); err != nil {
							fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
							continue
						}
						bufm_b[k] = complex(b, 0)
					} else if len(values) == 2 {
						var a, b float64
						if b, err = strconv.ParseFloat(values[0], 64); err != nil {
							fmt.Printf("String %s conversion to float error: %v\n", values[0], err)
							continue
						}
						if a, err = strconv.ParseFloat(values[1], 64); err != nil {
							fmt.Printf("String %s conversion to float error: %v\n", values[1], err)
							continue
						}
						bufm_a[k] = complex(a, 0)
						bufm_b[k] = complex(b, 0)
					}
					kk++
				}
				// Check for the normal EOF and the abnormal scan error
				// EOF does not give an error but is considered normal termination
				if !scanOK {
					if input.Err() != nil {
						fmt.Printf("Data file scan error: %v\n", input.Err().Error())
						return fmt.Errorf("data file scan error: %v", input.Err())
					}
				}
				// Put the next kk samples in back of the previous m
				if ARMA {
					copy(bufN_a[m:], bufm_a[:kk])
				}
				copy(bufN_b[m:], bufm_b[:kk])

				// window the m + kk samples with chosen window
				for i := 0; i < m+kk; i++ {
					if ARMA {
						bufN_a[i] *= w(i, m)
					}
					bufN_b[i] *= w(i, m)
				}

				// zero-pad N-m-kk samples in buf[N]
				for i := m + kk; i < N; i++ {
					if ARMA {
						bufN_a[i] = 0
					}
					bufN_b[i] = 0
				}

				// Perform N-point complex FFT and add squares to previous values in PSD[N/2+1]
				var fourierN_a []complex128
				if ARMA {
					fourierN_a = fft.FFT(bufN_a)
				}
				fourierN_b := fft.FFT(bufN_b)

				fb := cmplx.Abs(fourierN_b[0])
				if ARMA {
					fa := cmplx.Abs(fourierN_a[0])
					PSD[0] += fb * fb / (fa * fa)
				} else {
					PSD[0] += fb * fb
				}
				for i := 1; i < N/2; i++ {
					// Use positive and negative frequencies -> bufN[N-i] = bufN[-i]
					fb = cmplx.Abs(fourierN_b[i])
					fbNi := cmplx.Abs((fourierN_b[N-i]))
					if ARMA {
						fa := cmplx.Abs(fourierN_a[i])
						faNi := cmplx.Abs(fourierN_a[N-i])
						PSD[i] += (fa*fa + faNi*faNi) / (fb*fb + fbNi*fbNi)
					} else {
						PSD[i] += fb*fb + fbNi*fbNi
					}
				}
				fb = cmplx.Abs(fourierN_b[N/2])
				if ARMA {
					fa := cmplx.Abs(fourierN_a[N/2])
					PSD[N/2] += (fb * fb) / (fa * fa)
				} else {
					PSD[N/2] += fb * fb
				}

				// part of K*Sum(w[i]*w[i]) PSD normalizer
				normalizerPSD += sumWindow

				// EOF reached
				if !scanOK {
					break
				}
			} // K segments done

			// Normalize the PSD using K*Sum(w[i]*w[i])
			// Use log plot for wide dynamic range
			if r.FormValue("plottype") == "linear" {
				for i := range PSD {
					PSD[i] /= normalizerPSD
					if PSD[i] > psdMax {
						psdMax = PSD[i]
					}
					if PSD[i] < psdMin {
						psdMin = PSD[i]
					}
				}
				// log10 in dB
			} else {
				for i := range PSD {
					PSD[i] /= normalizerPSD
					PSD[i] = 10.0 * math.Log10(PSD[i])
					if PSD[i] > psdMax {
						psdMax = PSD[i]
					}
					if PSD[i] < psdMin {
						psdMin = PSD[i]
					}
				}
			}

			endpoints.xmin = 0
			endpoints.xmax = float64(N / 2) // equivalent to Nyquist critical frequency
			endpoints.ymin = psdMin
			endpoints.ymax = psdMax

			// Calculate scale factors for x and y
			xscale = (columns - 1) / (endpoints.xmax - endpoints.xmin)
			yscale = (rows - 1) / (endpoints.ymax - endpoints.ymin)

			// Store the PSD in the plot Grid
			for bin, pow := range PSD {
				row := int((endpoints.ymax-pow)*yscale + .5)
				col := int((float64(bin)-endpoints.xmin)*xscale + .5)
				plot.Grid[row*columns+col] = "online"
			}

			// Store in the form:  FFT Size, window type, number of samples nn, K segments, sampling frequency
			// Plot the PSD N/2 float64 values, execute the data on the plotfrequency.html template

			// Set plot status if no errors
			if len(plot.Status) == 0 {
				plot.Status = fmt.Sprintf("Status: Data plotted from (%.3f,%.3f) to (%.3f,%.3f)",
					endpoints.xmin, endpoints.ymin, endpoints.xmax, endpoints.ymax)
			}

		} else {
			// Set plot status
			fmt.Printf("Error opening file %s: %v\n", filename, err)
			return fmt.Errorf("error opening file %s: %v", filename, err)
		}
	} else {
		// Set plot status
		fmt.Printf("Error opening file %s: %v\n", filename, err)
		return fmt.Errorf("error opening file %s: %v", filename, err)
	}

	// Apply the  sampling rate in Hz to the x-axis using a scale factor
	// Convert the fft size to sampling rate/2, the Nyquist critical frequency
	sf := 0.5 * samplingRate / endpoints.xmax

	// Construct x-axis labels
	incr := (endpoints.xmax - endpoints.xmin) / (xlabels - 1)
	format := "%.0f"
	if incr*sf < 1.0 {
		format = "%.2f"
	}
	x := endpoints.xmin
	// First label is empty for alignment purposes
	for i := range plot.Xlabel {
		plot.Xlabel[i] = fmt.Sprintf(format, x*sf)
		x += incr
	}

	// Construct the y-axis labels
	incr = (endpoints.ymax - endpoints.ymin) / (ylabels - 1)
	y := endpoints.ymin
	for i := range plot.Ylabel {
		plot.Ylabel[i] = fmt.Sprintf("%.2f", y)
		y += incr
	}

	// Insert frequency domain parameters in the form
	plot.SampleFreq = fmt.Sprintf("%.0f", samplingRate)
	plot.FFTSegments = strconv.Itoa(K)
	plot.FFTSize = strconv.Itoa(N)
	plot.Samples = strconv.Itoa(nn)

	// Enter the filename in the form
	plot.Filename = path.Base(filename)

	// Write to HTTP using template and grid
	if err := freqTmpl.Execute(w, plot); err != nil {
		log.Fatalf("Write to HTTP output using template with error: %v\n", err)
	}

	return nil
}

// handlePlotResponse displays the impulse or frequency response of the
// LMS System Identification
func handlePlotResponse(w http.ResponseWriter, r *http.Request) {
	// main data structure
	var (
		plot PlotT
		err  error = nil
	)

	filename := r.FormValue("filename")
	// choose time or frequency domain processing
	if len(filename) > 0 {

		domain := r.FormValue("domain")
		switch domain {
		case "time":
			err = processTimeDomain(w, r, path.Join(dataDir, filename))
			if err != nil {
				plot.Status = err.Error()
			}
		case "frequency":
			err = processFrequencyDomain(w, r, path.Join(dataDir, filename))
			if err != nil {
				plot.Status = err.Error()
			}
		default:
			plot.Status = fmt.Sprintf("Invalid domain choice: %s", domain)
		}

		if err != nil {

			// Write to HTTP using template and grid``
			if err := plotresponseTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
		}
	} else {
		plot.Status = "Missing filename to plot"
		// Write to HTTP using template and grid``
		if err := plotresponseTmpl.Execute(w, plot); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}
	}
}

// generateSignal generates a new set of input data and sends it to the toChannel
func (lms *LMSAlgorithm) generateSignal() error {

	// generate a new set of input data for this trial

	// increment wg
	lms.wg.Add(1)

	// launch a goroutine to generate input samples consisting of
	// white Gaussian noise
	go func() {
		defer func() {
			lms.wg.Done()
			close(lms.toChan)
		}()

		// Send white Gaussian noise to channel via the toChannel
		for i := 0; i < lms.samples; i++ {
			// Send the sample through the channel
			lms.toChan <- rand.NormFloat64()
		}

	}()
	return nil
}

// generate the unknown Plant output and forward it along with the input noise to
// the next stage in the pipeline
func (lms *LMSAlgorithm) plant() error {
	lms.wg.Add(1)
	go func() {

		defer func() {
			lms.wg.Done()
			close(lms.fromChan)
		}()

		// number of coefficients in the numerator and denominator of the
		// plant difference equation, which is the order + 1
		ncoeff := 3 + 2*(len(lms.poleRad)-1)

		// Construct the biquads: 1-2*cos(zero-ang)z^(-1)+z^(-2)
		//                        -----------------------------------
		//                        1-2*r*cos(pole-ang)z^(-1)+r*r*z(-2)

		// each complex pole/zero pair is a biquad section
		nbiquads := len(lms.poleRad)
		biquad := make([]Biquad, nbiquads)

		// loop over the pole/zero pairs and construct the slice of biquads
		for i := 0; i < nbiquads; i++ {
			biquad[i] = Biquad{
				num:   [3]float64{1.0, -2.0 * math.Cos(lms.zeroAng[i]*deg2rad), 1.0},
				denom: [3]float64{1.0, -2.0 * lms.poleRad[i] * math.Cos(lms.poleAng[i]*deg2rad), lms.poleRad[i] * lms.poleRad[i]},
			}
		}

		// Compute the difference equation from the biquads by using polynomial convolution

		// resultant polynomials for numerator and denominator of the difference equation
		num := make([]float64, ncoeff)
		// initialize num with first biquad num
		for i := range num[0:3] {
			num[i] = biquad[0].num[i]
		}

		denom := make([]float64, ncoeff)
		// initialize denom with first biquad denom
		for i := range denom[0:3] {
			denom[i] = biquad[0].denom[i]
		}

		ncorr := 3 // number of correlations for this iteration in the for loop
		for i := 1; i < nbiquads; i++ {
			/***************** numerator ***************************************/
			num[0] += biquad[i].num[0] * num[0]
			num[1] += biquad[i].num[1]*num[0] + biquad[i].num[0]*num[1]
			for j := 2; j < ncorr; j++ {
				for k := 0; k < 3; k++ {
					num[j] += biquad[i].num[k] * num[j-k]
				}
			}
			num[ncorr] += biquad[i].num[2]*num[ncorr] + biquad[i].num[1]*num[ncorr+1]
			num[ncorr+1] += biquad[i].num[2] * num[ncorr+1]
			/***************** denominator ***************************************/
			denom[0] += biquad[i].denom[0] * denom[0]
			denom[1] += biquad[i].denom[1]*denom[0] + biquad[i].denom[0]*denom[1]
			for j := 2; j < ncorr; j++ {
				for k := 0; k < 3; k++ {
					denom[j] += biquad[i].denom[k] * denom[j-k]
				}
			}
			denom[ncorr] += biquad[i].denom[2]*denom[ncorr] + biquad[i].denom[1]*denom[ncorr+1]
			denom[ncorr+1] += biquad[i].denom[2] * denom[ncorr+1]

			ncorr += 2
		}

		// Save the plant model
		file, err := os.Create(path.Join(dataDir, plantModelFile))
		if err != nil {
			fmt.Printf("Create file %v error: %v", plantModelFile, err)
			close(lms.toChan)
			return
		}
		defer file.Close()
		for i := 0; i < ncoeff; i++ {
			fmt.Fprintf(file, "%f %f\n", num[i], denom[i])
		}

		// number of adaptive filter coefficients
		L := lms.order + 1
		// desired input to LMS algorithm is delayed by half the filter length
		delay := L / 2

		// x holds the previous inputs in the difference equation
		// y holds the previous ouputs in the difference equation
		x := make([]float64, L)
		y := make([]float64, L)

		// range over the to channel to obtain the white Gaussian noise from the generator
		i := 0
		for noise := range lms.toChan {
			out := 0.0
			// store the current input
			x[i] = noise
			for j := 0; j < ncoeff; j++ {
				k := (i - j) % L
				if k < 0 {
					k = L + k
				}
				out += num[j] * x[k]
			}

			for j := 1; j < ncoeff; j++ {
				k := (i - j) % L
				if k < 0 {
					k = L + k
				}
				out -= denom[j] * y[k]
			}

			// Send output to next stage via the from channel
			// "desired" is plant output delayed by L/2, "in" is the white noise
			// from the generator stage
			k := (i - delay) % L
			if k < 0 {
				k = L + k
			}
			lms.fromChan <- ComChan{desired: y[k], in: x[i]}

			// store the current output sample at the end of desired slice
			y[i] = out
			i = (i + 1) % L
		}
	}()
	return nil
}

// runLms determines the adaptive weights using the LMS algorithm
func (lms *LMSAlgorithm) runLms() error {

	// increment wg
	lms.wg.Add(1)

	// launch a goroutine to generate the adaptive filter
	go func() {
		defer lms.wg.Done()

		L := lms.order + 1
		// holds the previous inputs
		x := make([]float64, L)

		i := 0
		// range over the channel containing the plant output
		for d := range lms.fromChan {
			x[i] = d.in
			y := 0.0
			k := i
			for j := 0; j < L; j++ {
				y += lms.wTrial[j] * x[k]
				k = (k - 1) % L
				if k < 0 {
					k = L + k
				}
			}
			// calculate the error
			e := d.desired - y
			k = i
			// update the trial adaptive weights
			for j := 0; j < L; j++ {
				// normalize the update by the average square-
				// input multiplied by the number of filter coefficients
				lms.wTrial[j] += lms.gain * e * x[k] / float64(L)
				// update the ensemble weight
				lms.wEnsemble[j] += lms.wTrial[j]
				k = (k - 1) % L
				if k < 0 {
					k = L + k
				}
			}
			// Increment the current input index for x slice
			i = (i + 1) % L
		}
	}()

	return nil
}

// handleLmsSystemIdentifier creates a System Identifier using the LMS algorithm
func handleLmsSystemIdentifier(w http.ResponseWriter, r *http.Request) {
	var plot PlotT

	// Get the number of samples to generate and the sample rate
	// The number of samples determines the number of iterations of the LMS algorithm
	samplestxt := r.FormValue("samples")
	sampleratetxt := r.FormValue("samplerate")
	// choose time or frequency domain processing
	if len(samplestxt) > 0 && len(sampleratetxt) > 0 {

		samples, err := strconv.Atoi(samplestxt)
		if err != nil {
			plot.Status = fmt.Sprintf("Samples conversion to int error: %v", err.Error())
			fmt.Printf("Samples conversion to int error: %v", err.Error())
			// Write to HTTP using template and grid
			if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		samplerate, err := strconv.Atoi(sampleratetxt)
		if err != nil {
			plot.Status = fmt.Sprintf("Sample rate conversion to int error: %v", err.Error())
			fmt.Printf("Sample rate conversion to int error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		txt := r.FormValue("signalfreq")
		sigfreq, err := strconv.Atoi(txt)
		if err != nil {
			fmt.Printf("Signal frequency conversion error: %v\n", err)
			plot.Status = fmt.Sprintf("Signal frequency conversion to int error: %v", err.Error())
			// Write to HTTP using template and grid
			if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// Verify signal frequency is less than Nyquist frequency = sample rate / 2
		if sigfreq > samplerate/2 {
			fmt.Printf("Signal frequency %v is greater than Nyquist frequency %v\n", sigfreq, samplerate/2)
			plot.Status = fmt.Sprintf("Signal frequency %v is greater than Nyquist frequency %v\n", sigfreq, samplerate/2)
			// Write to HTTP using template and grid
			if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// adaptive filter order is the number of past samples
		// filter length = order + 1
		txt = r.FormValue("filtorder")
		order, err := strconv.Atoi(txt)
		if err != nil {
			plot.Status = fmt.Sprintf("Filter order conversion to int error: %v", err.Error())
			fmt.Printf("Filter order conversion to int error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}
		// Make filter length odd, which means filter order is even;
		// if order is 7, it is changed to 8, and the length becomes 9.
		order = (order + 1) / 2 * 2

		// gain factor that regulates the speed and stability of adaption
		txt = r.FormValue("gain")
		gain, err := strconv.ParseFloat(txt, 64)
		if err != nil {
			plot.Status = fmt.Sprintf("Gain conversion to float64 error: %v", err.Error())
			fmt.Printf("Gain conversion to float64 error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// number of trials to perform the LMS algorithm to get ensemble average of the
		// weights
		txt = r.FormValue("trials")
		trials, err := strconv.Atoi(txt)
		if err != nil {
			plot.Status = fmt.Sprintf("Trials conversion to int error: %v", err.Error())
			fmt.Printf("Trials conversion to int error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		txt = r.FormValue("snr")
		snr, err := strconv.Atoi(txt)
		if err != nil {
			plot.Status = fmt.Sprintf("SNR conversion to int error: %v", err.Error())
			fmt.Printf("SNR conversion to int error: %v\n", err.Error())
			// Write to HTTP using template and grid
			if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// Store the pole and zero data
		p := 1
		pr := make([]float64, 0)
		pa := make([]float64, 0)
		za := make([]float64, 0)
		for {
			txt = r.FormValue("polerad" + strconv.Itoa(p))
			if len(txt) == 0 {
				break
			}
			poleRad, err := strconv.ParseFloat(txt, 64)
			if err != nil {
				plot.Status = fmt.Sprintf("Pole radius conversion to float64 error: %v", err.Error())
				fmt.Printf("Pole radius conversion to float64 error: %v\n", err.Error())
				// Write to HTTP using template and grid
				if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}
			// Verify pole radius
			if poleRad < .01 || poleRad > .99 {
				plot.Status = fmt.Sprintf("Pole radius [%v] not between .01 and .99.\n", poleRad)
				fmt.Printf("Pole radius [%v] not between .01 and .99.\n", poleRad)
				// Write to HTTP using template and grid
				if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}
			pr = append(pr, poleRad)

			txt = r.FormValue("poleang" + strconv.Itoa(p))
			if len(txt) == 0 {
				break
			}
			poleAng, err := strconv.ParseFloat(txt, 64)
			if err != nil {
				plot.Status = fmt.Sprintf("Pole angle conversion to float64 error: %v", err.Error())
				fmt.Printf("Pole angle conversion to float64 error: %v\n", err.Error())
				// Write to HTTP using template and grid
				if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}

			// Verify pole angle
			if poleAng < 0 || poleAng > 180 {
				plot.Status = fmt.Sprintf("Pole angle [%v] not between 0 and 180.\n", poleAng)
				fmt.Printf("Pole angle [%v] not between 0 and 180.\n", poleAng)
				// Write to HTTP using template and grid
				if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}
			pa = append(pa, poleAng)

			txt = r.FormValue("zeroang" + strconv.Itoa(p))
			if len(txt) == 0 {
				break
			}
			zeroAng, err := strconv.ParseFloat(txt, 64)
			if err != nil {
				plot.Status = fmt.Sprintf("Zero angle conversion to float64 error: %v", err.Error())
				fmt.Printf("Zero angle conversion to float64 error: %v\n", err.Error())
				// Write to HTTP using template and grid
				if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}
			// Verify zero angle
			if zeroAng < 0 || zeroAng > 180 {
				plot.Status = fmt.Sprintf("Zero angle [%v] not between 0 and 180.\n", poleAng)
				fmt.Printf("Zero angle [%v] not between 0 and 180.\n", poleAng)
				// Write to HTTP using template and grid
				if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}
			za = append(za, zeroAng)
			p++
		}

		if len(pr) == 0 || len(pa) == 0 || len(za) == 0 {
			plot.Status = "Must have at least one pole and zero"
			fmt.Println("Must have at least one pole and zero")
			// Write to HTTP using template and grid
			if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}
		if (len(pr) != len(pa)) || (len(pr) != len(za)) || (len(pa) != len(za)) {
			plot.Status = "Must have the same number of pole radii, pole angles, and zero angles"
			fmt.Println("Must have the same number of pole radii, pole angles, an zero angles")
			// Write to HTTP using template and grid
			if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}

		// Construct object to hold LMS algorithm parameters
		lmsSystemIdentifier := LMSAlgorithm{
			gain:       gain,
			trials:     trials,
			order:      order,
			samples:    samples,
			samplerate: samplerate,
			wEnsemble:  make([]float64, order+1),
			wTrial:     make([]float64, order+1),
			snr:        snr,
			poleRad:    pr,
			poleAng:    pa,
			zeroAng:    za,
		}

		// Run the Least-Mean-Square (LMS) algorithm to create the adaptive filter
		// Loop over the trials to generate the ensemble of filters which is averaged.
		for i := 0; i < lmsSystemIdentifier.trials; i++ {
			// Create new channels each trial since they are closed after each trial
			lmsSystemIdentifier.fromChan = make(chan ComChan)
			lmsSystemIdentifier.toChan = make(chan float64)

			// Generate a new set of input data consisting of a sine wave at unit amplitude
			err = lmsSystemIdentifier.generateSignal()
			if err != nil {
				plot.Status = fmt.Sprintf("generateSignal error: %v", err.Error())
				fmt.Printf("generateSignal error: %v\n", err.Error())
				// Write to HTTP using template and grid
				if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}

			// Run the signal through the unknown Plant
			err = lmsSystemIdentifier.plant()
			if err != nil {
				plot.Status = fmt.Sprintf("chanFilter error: %v", err.Error())
				fmt.Printf("chanFilter error: %v\n", err.Error())
				// Write to HTTP using template and grid
				if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}

			// Run the LMS algorithm to find the adaptive filter coefficients
			err = lmsSystemIdentifier.runLms()
			if err != nil {
				plot.Status = fmt.Sprintf("generateLMSData error: %v", err.Error())
				fmt.Printf("run LmsSystemIdentifer error: %v\n", err.Error())
				// Write to HTTP using template and grid
				if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
					log.Fatalf("Write to HTTP output using template with error: %v\n", err)
				}
				return
			}
			// Wait for this trial to finish
			lmsSystemIdentifier.wg.Wait()
		}

		// Save ensemble filter coefficients to disk (averaged coefficients)
		file, err := os.Create(path.Join(dataDir, lmsSystemIdentifierFile))
		if err != nil {
			plot.Status = fmt.Sprintf("create %v error: %v", path.Join(dataDir, lmsSystemIdentifierFile), err.Error())
			fmt.Printf("create %v error: %v\n", path.Join(dataDir, lmsSystemIdentifierFile), err.Error())
			// Write to HTTP using template and grid
			if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
				log.Fatalf("Write to HTTP output using template with error: %v\n", err)
			}
			return
		}
		defer file.Close()

		// Average the weight ensemble to reduce variance
		for i := 0; i < lmsSystemIdentifier.order+1; i++ {
			fmt.Fprintf(file, "%v\n", lmsSystemIdentifier.wEnsemble[i]/(float64(trials*samples)))
		}

		plot.Status = fmt.Sprintf("Adaptive filter weights '%s' written to the data folder", lmsSystemIdentifierFile)
		// Write to HTTP using template and grid
		if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}

	} else {
		plot.Status = "Enter samples, sample rate, SNR, etc."
		// Write to HTTP using template and grid
		if err := lmsSystemIdentifierTmpl.Execute(w, plot); err != nil {
			log.Fatalf("Write to HTTP output using template with error: %v\n", err)
		}
	}
}

// executive program
func main() {

	// Setup http server with handler for creating a System Identifier using LMS
	http.HandleFunc(patternLmsSystemIdentifier, handleLmsSystemIdentifier)

	// Setup http server with handler for plotting impulse or frequency responses for the LMS System Identifier
	http.HandleFunc(patternPlotResponse, handlePlotResponse)

	fmt.Printf("LMS System Identifier Server listening on %v.\n", addr)

	http.ListenAndServe(addr, nil)
}
