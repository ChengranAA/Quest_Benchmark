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

![Screenshot 2025-02-18 at 1.56.03 PM](/Users/lcraaaa/Desktop/quest_extension/demo/README.assets/Screenshot 2025-02-18 at 1.56.03 PM.png)

## API

Some C++ APIs are hidden from the wrapper but it is nice to have for debugging. 

- `get_pdf` returns the current probablity density function into an Godot Array 
  - `print(quest._quest.get_pdf())`
- `get_intensity_history` returns the intensity history into an Godot Array 
  - `print(quest._quest.get_intensity_history())`
- `get_response_history` returns the response history into an Godot Array 
  - `print(quest._quest.get_reponse_history())`

## Bug 

Record this into the GitHub issue 
