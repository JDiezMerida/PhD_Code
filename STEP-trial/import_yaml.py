import yaml
import os
import copy

folder= 'C:/Users/jdiez/Documents/GitHub/PhD_Code/STEP-trial/docs/'
# Load all YAML files
def load_protocols(folder=folder):
    protocols = {}
    for file in os.listdir(folder):
        if file.endswith(".yaml"):
            with open(os.path.join(folder, file), "r") as f:
                data = yaml.safe_load(f)
                pid = data["protocol"]["id"]
                protocols[pid] = data
    return protocols


# Merge child with parent
def resolve_inheritance(protocols, proto):
    if "inherits_from" not in proto:
        return proto

    parent_id = proto["inherits_from"]
    parent = protocols[parent_id]

    # Recursively resolve parent
    parent_resolved = resolve_inheritance(protocols, parent)

    # Deep copy parent
    merged = copy.deepcopy(parent_resolved)

    # Apply modifications
    if "modifications" in proto:
        for mod in proto["modifications"]:
            step_id = mod["step_id"]

            for section in merged["sections"]:
                for step in section["steps"]:
                    if step["id"] == step_id:
                        step.update(mod)

    # Override metadata if present
    if "metadata" in proto:
        merged["metadata"] = proto["metadata"]

    # Override title if present
    if "protocol" in proto:
        merged["protocol"].update(proto["protocol"])

    return merged


# Render to Markdown
def render(protocol):
    out = []

    out.append(f"# {protocol['protocol']['title']}\n")

    meta = protocol.get("metadata", {})
    if meta:
        out.append(f"**Authors:** {', '.join(meta.get('authors', []))}")
        out.append(f"**Created:** {meta.get('created', '')}")
        out.append(f"**Last modified:** {meta.get('last_modified', '')}\n")

    if "description" in protocol:
        out.append("## Description\n")
        out.append(protocol["description"] + "\n")

    for section in protocol["sections"]:
        out.append(f"# {section['title']}\n")

        out.append("| Step | Task | Subtask | Materials | Parameters | Warnings | Verification |")
        out.append("|------|------|-------|------------|----------|--------------|---------|")

        for step in section["steps"]:
            row = f"| {step.get('id','')} | {step.get('task','')} | {step.get('subtask','')} | {step.get('materials','')} | {step.get('parameters','')} | {step.get('warnings','')} | {step.get('verification','')} |"
            out.append(row)

        out.append("\n")

    return "\n".join(out)


# Main
protocols = load_protocols()


OUTPUT_DIR = folder

for pid, proto in protocols.items():
    resolved = resolve_inheritance(protocols, proto)
    content = render(resolved)

    # Save inside folder based on category (default: root)
    folder = resolved.get("category", "general")

    folder_path = os.path.join(OUTPUT_DIR, folder)
    os.makedirs(folder_path, exist_ok=True)

    filename = f"{pid}.md"
    filepath = os.path.join(folder_path, filename)

    with open(filepath, "w") as f:
        f.write(content)

'''
# Choose which protocol to render
TARGET = "exfoliation_hbn_v1"

resolved = resolve_inheritance(protocols, protocols[TARGET])

dir_list= os.listdir(folder)
#print(dir_list)
filename = [f for f in dir_list if f.endswith(".yaml")]

for i in filename:
    with open(folder+i[:-5]+".md", "w") as f:
        f.write(render(resolved))
'''