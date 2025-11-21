/**
 * Page Transition Component - Smooth page transitions with animations
 */
import React, { useEffect, useState } from 'react'
import { useLocation } from 'react-router-dom'

interface PageTransitionProps {
  children: React.ReactNode
  duration?: number
}

export const PageTransition: React.FC<PageTransitionProps> = ({
  children,
  duration = 300,
}) => {
  const location = useLocation()
  const [displayLocation, setDisplayLocation] = useState(location)
  const [transitionStage, setTransitionStage] = useState<'fadeIn' | 'fadeOut'>('fadeIn')

  useEffect(() => {
    if (location !== displayLocation) {
      setTransitionStage('fadeOut')
    }
  }, [location, displayLocation])

  useEffect(() => {
    if (transitionStage === 'fadeOut') {
      const timeout = setTimeout(() => {
        setDisplayLocation(location)
        setTransitionStage('fadeIn')
      }, duration)

      return () => clearTimeout(timeout)
    }
  }, [transitionStage, location, duration])

  return (
    <div
      className={`page-transition ${transitionStage}`}
      style={{
        animation: `${transitionStage} ${duration}ms ease-in-out`,
      }}
    >
      {children}
    </div>
  )
}

export default PageTransition
