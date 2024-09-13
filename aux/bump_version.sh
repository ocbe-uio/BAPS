#!/bin/bash

# Read the file line by line
while IFS= read -r line
do
  # Check if the line contains the version vector
  if [[ $line == *"ver = "* ]]; then
    # Extract the version vector
    version=$(echo $line | cut -d'=' -f2 | tr -d '[];')

    # Split the version vector into an array
    IFS=' ' read -r -a versionArray <<< "$version"

    # Increment the last element of the version vector
    ((versionArray[${#versionArray[@]}-1]++))

    # Reconstruct the version vector
    version=$(IFS=' '; echo "${versionArray[*]}")

    # Reconstruct the line
    line="  ver = [$version];"
  fi

  # Print the line
  echo "$line"
done < "baps.m" > "baps_new.m"


# Replace the old file with the new file
mv "baps_new.m" "baps.m"

# Final message (useful for commit)
version_string=$(echo "$version" | tr ' ' '.')
echo "Updated BAPS version to $version_string"
