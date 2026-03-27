import yaml
folder= 'C:/Users/jdiez/Documents/GitHub/PhD_Code/STEP-trial/docs/'
with open(folder+"Graphene_exfoliation.yaml", "r") as f:
    data = yaml.safe_load(f)

out = []

# Title
out.append(f"# {data['protocol']['title']}\n")

# Metadata
meta = data["metadata"]
out.append(f"**Authors:** {', '.join(meta['authors'])}")
out.append(f"**Created:** {meta['created']}")
out.append(f"**Last modified:** {meta['last_modified']}\n")

# Description
out.append("## Description\n")
out.append(data["description"] + "\n")

# Sections
for section in data["sections"]:
    out.append(f"# {section['title']}\n")

    out.append("| Step | Name | Materials | Parameters | Warnings | Verification |")
    out.append("|------|------|-----------|------------|----------|--------------|")

    for step in section["steps"]:
        row = f"| {step['id']} | {step['name']} | {step['materials']} | {step['parameters']} | {step['warnings']} | {step['verification']} |"
        out.append(row)

    out.append("\n")

# Save markdown
with open(folder+"Graphene_exfoliation.md", "w") as f:
    f.write("\n".join(out))