<h4>LMS System Identifier and Modeler</h4>

<p>
This program is a web application written in Go which uses the html/template package to create the HTML.  To start the server, issue "go run identification.go" in the Command Prompt.  To create an adaptive system identifier/modeler using LMS, enter http://127.0.0.1:8080/lmssystemidentifier in the web browser address bar.  Enter the <i>LMS System Identifier Options</i> form data.  The adaptive filter coefficients are averaged over the specified number of trials to produce an ensemble average in order to reduce the variance. Enter the number of samples, sample rate, LMS gain, trials, and adaptive filter order.  The gain is a small number about 0.001 or less.  The pole and zero table specifies the unknown plant system which the adaptive filter will identify.  Not all entries will produce a stable plant.  If the impulse response of the plant is increasing over time it means the plant is unstable and cannot be modeled by the adaptive filter.  Another set of poles and zeros will need to be tried in order to create a stable plant that does not "blow up."  The poles and zeroes are located inside the unit circle in the complex z-plane.  Each pole and zero is a complex conjugate pair.  Up to five complex conjugate poles or zeros can be specified for the plant.
</p>
<p>
The poles and zeros are converted into quadratic pole/zero stages, or <i>biquads</i>.  The biquads are multiplied together using polynomial multiplication, which is a convolution sum.  The resultant polynomial is the system or transfer function or z-transform of the plant.  The difference equation can be determined by inspection of the transfer function.  The difference equation is used to compute the impulse response of the plant by supplying a unit sample as the input to the difference equation.  The adaptive filter contains no poles; it is a FIR filter.  The source to both the plant and the adaptive filter is white Gaussian noise with unit variance.  The output of the plant going to the LMS algorithm error detector is delayed by half the length of the adaptive filter in order to make the adaptive filter causal and have linear phase response.
</p>
<p>
After the adaptive filter is produced, the impulse or frequency response of the plant or adaptive filter can be plotted.  The link <i>Plot Response</i> will take you to the page where you can choose the impulse (time domain) or frequency response.  For the frequency response, the FFT size, number of segments, and window can be chosen.  For this program, one segment with a rectangular window and an FFT size of 8192 is sufficient.
</p>

<h4>LMS System Identifier Options, Bandpass, Order 80 Adaptive Filter</h4>

![LmsSystemIdentifierOptions](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/5140f3ab-92d0-493d-9ef3-a28691ffc149)

<h4>Plant Impulse Response</h4>

![impulseResponsePlant](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/5043f982-5921-4926-973e-cbb4f43b10a2)

<h4>Plant Frequency Response</h4>

![freqResponsePlant](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/2bff1f93-1131-4186-8479-5fb29da155b5)

<h4>Adaptive Filter Impulse Response, Order 80</h4>

![impulseResponseAdaptiveFilt](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/91abf45e-519b-4ca1-a1b8-899fdbbfa491)

<h4>Adaptive Filter Frequency Response, Order 80</h4>

![freqResponseAdaptiveFilt](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/a05a8ce0-39ba-48a4-835f-2d6229431629)

<h4>LMS System Identifier Options, Lowpass, Order 80 Adaptive Filter</h4>

![LmsSystemIdentifierOptions](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/d40ca8c4-29ac-4297-a6be-5da61248799c)

<h4>Plant Impulse Response</h4>

![impulseResponsePlant](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/79940224-c4ae-48ad-a4ec-d6735d4816f2)

<h4>Plant Frequency Response</h4>

![frequencyResponsePlant](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/512efbd7-60be-485b-8eb2-3ece4b51ade6)

<h4>Adaptive Filter Impulse Response, Order 80</h4>

![impulseResponseAdaptiveFilter](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/43cc65a4-db13-4d77-8911-91659a8d0fe7)

<h4>Adaptive Filter Frequency Response, Order 80</h4>

![frequencyResponseAdaptiveFilter](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/b236226b-5272-4c0c-950e-6ee303b00d0a)

<h4>LMS System Identifier Options, Highpass, Order 80 Adaptive Filter</h4>

![LmsSystemIdentifierOptions](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/3b67b487-c183-44f1-be89-9cfe387d44f9)

<h4>Plant Impulse Response</h4>

![impulseResponsePlant](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/5100566b-fd7c-4f2e-b082-5816b3d6d62c)

<h4>Plant Frequency Response</h4>

![frequencyResponsePlant](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/790d48ec-dc35-40c2-a8e4-34c5da4ba7f7)

<h4>Adaptive Filter Impulse Response, Order 80</h4>

![impulseResponseAdaptiveFilter](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/0670fde8-2eb1-40f6-bb71-6f3040993de9)

<h4>Adaptive Filter Frequency Response, Order 80</h4>

![frequencyResponseAdaptiveFilter](https://github.com/thomasteplick/LmsSystemIdentification/assets/117768679/db9881e5-46a6-4522-836e-3cd56d2f9b2d)












