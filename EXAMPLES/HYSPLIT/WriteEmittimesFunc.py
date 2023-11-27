def find_line_numbers(file_path, target_string):
    line_numbers=[]
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        total_lines=len(lines) + 1
        # Iterate over the lines and write line numbers for lines containing the target string
        for line_number, line_content in enumerate(lines, start=1):
            if target_string in line_content:
                line_numbers.append(line_number)
    except FileNotFoundError:
        print(f"Error: File not found - {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")
    line_numbers.append(total_lines)
    return(line_numbers)

def compute_number_of_lines(integer_list):
    differences = []

    for i in range(len(integer_list) - 1):
        difference = (integer_list[i + 1] - integer_list[i]) -1
        differences.append(difference)

    return differences


def process_file(file_path, target_string, replacement_list):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()

        count_occurrences = 0

        # Iterate over the lines
        for i, line in enumerate(lines):
            if target_string in line:
                count_occurrences += 1
                # Replace the occurrence with the corresponding element from the list
                replacement_index = count_occurrences - 1
                if replacement_index < len(replacement_list):
                    replacement_value = replacement_list[replacement_index]
                    lines[i] = line.replace(target_string, str(replacement_value))

        # Write the modified lines back to the file
        with open(file_path, 'w') as file:
            file.writelines(lines)

    except FileNotFoundError:
        print(f"Error: File not found - {file_path}")
    except Exception as e:
        print(f"An error occurred: {e}")


def write_emittimes(file_path, output_file_path, target_element_count,max_rows,vent_lat,vent_lon,vent_height):
        
         
        output_file=open(output_file_path,'a')
        file=open(file_path,'r')
         
        for line in file:
        
            # Split the line into elements
            elements = line.strip().split()
            # Check if the number of elements matches the target
            if len(elements) == target_element_count:
                # If a row with the target element count is found, start collecting
                current_row = elements
                num_following_rows = int(elements[5])
                current_row[5] = str(max_rows)
                following_rows = []

                for _ in range(num_following_rows):
                    following_line = file.readline()
                    if following_line:
                        following_elements = following_line.strip().split()
                        following_rows.append(following_elements)
                        
                if len(following_rows) < max_rows:
                    rowtoappend = following_elements                
                    rowtoappend[6] = str(vent_lat)
                    rowtoappend[7] = str(vent_lon)
                    rowtoappend[8] = str(vent_height)
                    rowtoappend[9] = str(0.000)
                    rowtoappend[10] = str(0.000)
                    for index in range(int(max_rows-len(following_rows))):
                        following_rows.append(rowtoappend)
                               
                # Write current_row and following_rows to the output file
                output_file.write("{}\n".format(" ".join(current_row)))
                for row in following_rows:
                    output_file.write("{}\n".format(" ".join(row)))
                            
