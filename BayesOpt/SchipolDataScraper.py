from bs4 import BeautifulSoup
from urllib.request import urlopen
import time
import matplotlib.pyplot as plt


# Taken from https://www.codegrepper.com/code-examples/python/days+in+month+function+python
def leap_year(year):
    if year % 400 == 0:
        return True
    if year % 100 == 0:
        return False
    if year % 4 == 0:
        return True
    return False


def days_in_month(month, year):
    if month in {1, 3, 5, 7, 8, 10, 12}:
        return 31
    if month == 2:
        if leap_year(year):
            return 29
        return 28
    return 30


start = time.time()
AircraftList = []
for month in range(1, 13):
    monthStr = str(month)
    if len(monthStr) == 1:
        monthStr = '0' + monthStr
    dayNumber = days_in_month(month, 2018)
    for day in range(1, dayNumber + 1):
        dayStr = str(day)
        if len(dayStr) == 1:
            dayStr = '0' + dayStr
        url = "https://schiphol.dutchplanespotters.nl/?date=2018-" + monthStr + "-" + dayStr + "&group=aircraft"
        print(url)

        page = urlopen(url)
        html = page.read().decode("utf-8")
        soup = BeautifulSoup(html, "html.parser")

        # Generate Table with Information
        searchString = ["day0", "day1", "gold0", "gold1", "dark0 noprint", "dark1 noprint", "day0wishlist",
                        "day1 wishlist", "dark0 wishlist", "dark1 wishlist", "dark0 noprint wishlist",
                        "dark1 noprint wishlist"]
        for string in searchString:
            table = soup.body.table.tbody.findAll("tr", {"class": string})
            for subTable in table:
                subTable = subTable.get_text()
                subTable = subTable.split("\n")
                subTable.remove('')
                subTable.remove('')
                if subTable[0] == 'Cargo':
                    AircraftList.append([subTable[1], subTable[2]])
                else:
                    AircraftList.append([subTable[0], subTable[1]])

halfway = time.time()
print(halfway - start)
aircraftEvents = []
totalCount = 0
while True:
    AircraftNumberList = []
    aircraftCombination = AircraftList[0]
    count_a = AircraftList.count(aircraftCombination)
    totalCount += count_a
    aircraftEvents.append([count_a, aircraftCombination[0], aircraftCombination[1]])
    AircraftList = [value for value in AircraftList if value != aircraftCombination]
    if len(AircraftList) == 0:
        break
aircraftEvents = sorted(aircraftEvents, key=lambda x: x[0], reverse=True)

print(aircraftEvents)
partialCount = 0
percentRule = 0.98 * totalCount
reqMonths = []
cumPerc = []
for element in aircraftEvents:
    print(element)
    partialCount += element[0]
    monthsRequired = 458.3 / (element[0] * 1.5)
    print("months required", monthsRequired, "cumulative percentage flights", partialCount / totalCount * 100)
    reqMonths.append(monthsRequired)
    cumPerc.append(partialCount / totalCount * 100)
    if partialCount > percentRule or monthsRequired > 60:
        break
print(totalCount)
end = (halfway - time.time())

plt.plot(reqMonths, cumPerc)
plt.title("Estimated Percentage of Aircraft Type-Airline Combinations Modelled against Time")
plt.xlabel("Time [months]")
plt.ylabel("Percentage Modelled [%]")
plt.show()

