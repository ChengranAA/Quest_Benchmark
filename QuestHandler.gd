class_name QuestHandler

# things should inherited from parent class StairHandler 
var data: Array
var nTrials
var thisTrialN = -1
var finished: bool = false

var startVal        # Initial guess for stimulus contrast (mid-range)
var startValSd : float          # Standard deviation of prior (initial uncertainty)
var pThreshold : float           # Target performance level
var method 
var beta                      # Slope of the psychometric function
var delta                    # Lapse rate
var gamma                    # Guess rate
var grain
var minVal
var maxVal
var StopInterval
var _questNextIntensity
var _nextIntensity
var _range 

var _quest = Quest.new()
var intensities: Array

func quantile():
	return _quest.quantile()
	

func _intensity():
	self._questNextIntensity = self.quantile()
	self._nextIntensity = self._questNextIntensity
	

# this is called with the QuestHandler.new(args) 
func _init(startVal, startValSd, pThreshold,nTrials, 
			stopInterval = null, method = 'quantile', 
			beta = 3.5, delta = 0.01, gamma = 0.5, grain = 0.01, 
			range = null, minVal = null, maxVal = null):
				
	self.startVal = startVal
	self.startValSd = startValSd
	self.pThreshold = pThreshold
	self.nTrials = nTrials
	self.StopInterval = stopInterval # not implemneted don't use 
	self.method = method # only quantile supported by the quest object  
	self.beta = beta
	self.delta = delta
	self.gamma = gamma
	self.grain = grain
	self._range = range
	self.minVal = minVal
	self.maxVal = maxVal
	
	self._questNextIntensity = startVal
	
	_quest.initialize(startVal, startValSd, pThreshold, beta, delta, gamma)
	
func next():
	if self.finished == false: 
		self.thisTrialN += 1
		self.intensities.append(self._nextIntensity)
		return self._nextIntensity
	else:
		self._terminate()


func _checkFinished():
	if self.nTrials != null and len(self.intensities) >= self.nTrials:
		self.finished = true
	else: 
		self.finished = false

	
func calculateNextIntensity():
	# based on current intensity and counter of correct responses
	self._intensity()
	if self.maxVal != null and self._nextIntensity > self.maxVal: 
		self._nextIntensity = self.maxVal
	elif self.minVal != null and self._nextIntensity < self.minVal:
		self._nextIntensity = self.minVal
	self._questNextIntensity = self._nextIntensity


func addResponse(result, intensity = null):
	
	# Add a 1 or 0 to signify a correct / detected or
	# incorrect / missed trial
	# Supplying an `intensity` value here indicates that you did not use the
	# recommended intensity in your last trial and the staircase will
	# replace its recorded value with the one you supplied here.
		
	if intensity == null:
		intensity = self._questNextIntensity
	else: 
		if len(self.intensities) != 0:
			self.intensities.pop_back()
		self.intensities.push_back(intensity) 
	
	# update quest 
	self._quest.update(intensity, result)
	self.data.push_back(result)
	
	self._checkFinished() 
	if not self.finished:
		self.calculateNextIntensity()
		
func _terminate():
	pass
	

	
