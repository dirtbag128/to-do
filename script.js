document.addEventListener("DOMContentLoaded", function () {
  const body = document.querySelector("body");
  let gradientIndex = 0;
  const gradients = [
    "linear-gradient(45deg, #ff7e5f, #feb47b)", // Soft orange-pink gradient
    "linear-gradient(45deg, #6a11cb, #2575fc)", // Purple-blue gradient
    "linear-gradient(45deg, #ff9a8b, #ffc3a0)", // Peach-pink gradient
    "linear-gradient(45deg, #24c6dc, #514a9d)", // Blue-purple gradient
  ];

  function changeBackground() {
    body.style.background = gradients[gradientIndex];
    gradientIndex = (gradientIndex + 1) % gradients.length; // Loop through gradients
  }
  setInterval(changeBackground, 5000); // Change background every 5 seconds
});
