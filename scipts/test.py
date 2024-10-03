import os

# Create the logs directory if it doesn't exist
os.makedirs('/beegfs/scratch/ric.broccoli/kubacki.michal/logs', exist_ok=True)

# Write the message to a text file in the logs directory
with open('/beegfs/scratch/ric.broccoli/kubacki.michal/logs/hello_world.txt', 'w') as f:
    f.write("Hello World")

print("Message saved to /beegfs/scratch/ric.broccoli/kubacki.michal/logs/hello_world.txt")