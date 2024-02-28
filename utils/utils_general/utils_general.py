def read_text_file(file_path):
    lines = []
    with open(file_path, 'r') as file:
        for line in file:
            cleaned_line = line.strip()
            if not cleaned_line.startswith("#"):
                lines.append(cleaned_line)
    return lines


