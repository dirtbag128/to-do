document.addEventListener("DOMContentLoaded", function () {
  // 1. Simple gradient background animation
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

  // 2. Animate headings and tasks with subtle effects
  const headings = document.querySelectorAll(".todo-list h1, .todo-list h2");
  const tasks = document.querySelectorAll(".task-list li");

  // Animate headings with fade-in and slight upward movement
  headings.forEach((heading, index) => {
    heading.style.opacity = 0;
    heading.style.transform = "translateY(-20px)";
    heading.style.transition = "opacity 1s ease, transform 1s ease";
    setTimeout(() => {
      heading.style.opacity = 1;
      heading.style.transform = "translateY(0)";
    }, index * 300); // Delay based on index for staggered effect
  });

  // Animate tasks with fade-in and slight leftward movement
  tasks.forEach((task, index) => {
    task.style.opacity = 0;
    task.style.transform = "translateX(-20px)";
    task.style.transition = "opacity 0.8s ease, transform 0.8s ease";
    setTimeout(
      () => {
        task.style.opacity = 1;
        task.style.transform = "translateX(0)";
      },
      headings.length * 300 + index * 100,
    ); // Staggered delay for tasks
  });
});
