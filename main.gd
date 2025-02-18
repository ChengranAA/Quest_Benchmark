extends Node3D

const QuestHandlerMod = preload("res://QuestHandler.gd") # this load the quest handler class 

# Quest Setting 
var start_contrast = 0.5 
var contrast_sd = 0.3
var pThreshold = 0.5
var beta = 3.5 
var delta = 0.01
var gamma = 0.01

var nTrials = 30
var contrast_range = [0.01, 1.0]

var quest = QuestHandlerMod.new(start_contrast, 
								contrast_sd, 
								pThreshold, 
								nTrials, 
								null, 
								"quantile", 
								beta, 
								delta,
								gamma, 
								0.01, 
								null,
								contrast_range[0], 
								contrast_range[1])


func _ready():
	var new_intensity
	quest.addResponse(1)
	new_intensity = quest.next()
	print(new_intensity)
	print(quest._quest.get_pdf())
	
	

# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta: float) -> void:
	pass
