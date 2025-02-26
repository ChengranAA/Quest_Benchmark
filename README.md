# Example Project for Staircasing 

## Project Steup

Install godot 4.3 (stable)

Either from the official website https://godotengine.org or from brew 

```zsh
brew install --cask godot 
```

Clone this project 

```zsh
git clone <this project url>
```

Open the app and import the project 

Open the little paper icon to adjust the script, if you don't see the editor open, then click the script tab. 

![Screenshot 2025-02-18 at 1.56.03 PM](https://github.com/ChengranAA/Quest_Benchmark/blob/main/README.assets/Screenshot%202025-02-18%20at%201.56.03%E2%80%AFPM.png)

Scipts for testing should be in the `main.gd` file and under the `func _ready()`. 

## API

Some C++ APIs are hidden from the wrapper but it is nice to have for debugging. 

- `get_pdf` returns the current probablity density function into an Godot Array 
  - `print(quest._quest.get_pdf())`
- `get_intensity_history` returns the intensity history into an Godot Array 
  - `print(quest._quest.get_intensity_history())`
- `get_response_history` returns the response history into an Godot Array 
  - `print(quest._quest.get_reponse_history())`

## Bug 

Record bugs into the GitHub issue 
