# SPDX-License-Identifier: GPL-3.0-or-later
# Droplet Detective
# Project: https://github.com/3DCP-TUe/DropletDetective
#
# Copyright (c) 2024 Eindhoven University of Technology
#
# Authors:
#   - Arjen Deetman (2024)
#
# For license details, see the LICENSE file in the project root.

import asyncio
import logging
import csv
import os
from datetime import datetime
from asyncua import Client, ua

async def main():
  
    while True:

        client = Client(url="opc.tcp://10.129.4.2:4840")
        client.set_user("Admin") 
        client.set_password("admin")

        database = "D:/GitHub/DropletDetective/20240528_Test.csv"
        create_database(database)

        counter = 0

        try:
            async with client:

                node = client.get_node("ns=1;i=104") 
                
                while True: 
                    
                    # Time
                    date = datetime.now()
                    t = date.strftime("%H:%M:%S.%f")[:-3]

                    # Value
                    value = await node.read_value()

                    # Write data to CSV database
                    with open(database, 'a', newline='') as file:
                        writer = csv.writer(file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                        writer.writerow([t, "{0:.8f}".format(value)])
                    
                    if counter > 100:
                        logging.info("{0}, {1:.8f}".format(t, value))
                        counter = 0
                    
                    counter += 1

                    await asyncio.sleep(0.020)
                    await client.check_connection() # Throws a exception if connection is lost

        except ua.UaError as e:
            logging.warning("An OPC UA error occurred: {}".format(e))
        except ConnectionError:
            logging.warning("Lost connection to the OPC UA server. Reconnecting in 2 seconds...")
            await asyncio.sleep(2)


def create_database(file: str) -> None:
        
    """Sets and initiates the CSV database."""

    # Initiate file with header (only if file doesn't exist yet)
    if not (os.path.exists(file)):
        with open(file, 'a', newline='') as file:
            writer = csv.writer(file, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["Time", "Load"])


if __name__ == "__main__":
    logging.getLogger('asyncua').setLevel(logging.WARNING)
    logging.basicConfig(level=logging.INFO)
    asyncio.run(main())
